suppressPackageStartupMessages(library(package = "RCurl"))
suppressPackageStartupMessages(library(package = "readr"))
suppressPackageStartupMessages(library(package = "gdata"))
suppressPackageStartupMessages(library(package = "xlsx"))
suppressPackageStartupMessages(library(package = "affy"))
suppressPackageStartupMessages(library(package = "hgu133a2cdf"))
suppressPackageStartupMessages(library(package = "plyr"))
suppressPackageStartupMessages(library(package = "dplyr"))
suppressPackageStartupMessages(library(package = "ggplot2"))
suppressPackageStartupMessages(library(package = "tidyr"))
suppressPackageStartupMessages(library(package = "tibble"))
suppressPackageStartupMessages(library(package = "synapseClient"))
suppressPackageStartupMessages(library(package = "jsonlite"))
suppressPackageStartupMessages(library(package = "httr"))
suppressPackageStartupMessages(library(package = "XML"))
suppressPackageStartupMessages(library(package = "igraph"))
suppressPackageStartupMessages(library(package = "org.Hs.eg.db"))
suppressPackageStartupMessages(library(package = "SuperExactTest"))


# set default options/variables
options(stringsAsFactors = FALSE)
options(useFancyQuotes = FALSE)

# convert probeIDs to gene symbol
# remove probes matching multiple genes
# for many probes matching one gene - keep the probe that has the maximum median across all samples
source("probe2gene_V2.R")

# read annotation and expression data
hgu133a_2 <- data.frame(read_csv("HG-U133A_2.na35.annot.csv", comment = "#"))
expression_data <- read_tsv("ViralChallenge_training_EXPRESSION_RMA.tsv", 
												comment = "#") %>%
							   as.data.frame()
# fix the dimention and probe names
row.names(expression_data) <- expression_data$FEATUREID
expression_data <- expression_data %>% dplyr::select(-FEATUREID)
expression_data <- t(expression_data) %>% as.data.frame()
    
# fix the column name
colnames(hgu133a_2)[1] <- "ID"

# save
expression_data_gene_level <- megre_probe2gene(GE_df = expression_data, 
																				    GPL = hgu133a_2)

# make annotation data frame
annotDF <- data.frame(colNames = colnames(expression_data_gene_level)) %>%
				   mutate(colNames = gsub(pattern = "^X",
				   											 replacement = "",
				   											 colNames), 
				   				PROBEID = gsub(pattern = "(+):.+",
				   										     replacement = "\\1",
				   										     colNames),
				   			   GENE = gsub(pattern = ".+:(.+)",
				   									   replacement = "\\1",
				   									   colNames)) %>%
				   				dplyr::select(-colNames)
rownames(annotDF) <- annotDF$PROBEID			

# clinical data
clinical_data <- read_tsv("ViralChallenge_training_CLINICAL.tsv", comment = "#") %>%
						 as.data.frame() %>%
						 mutate(TIME = ifelse(TIMEHOURS <=0, "TIME0", "TIME24"),
						 			 VIRUS = ifelse(STUDYID %in% STUDYID[grep("H1N1", STUDYID)],
						 			 						  "H1N1",
						 			 				ifelse(STUDYID %in% STUDYID[grep("H3N2", STUDYID)],
						 			 						  "H3N2", "Rhinovirus")))

# calculate fold change of genes (1 vs 0 for SC1 and SC2, and ~SC3)
# for SC1 and SC2, calculate log2(FC) : 1 vs 0
matDF <- expression_data_gene_level
colnames(matDF) <- gsub(pattern = ".+:(.+)", 
										   replacement = "\\1",
										   colnames(matDF))
matDF <- matDF %>%
				rownames_to_column()
matDF1 <- matDF %>%
				  mutate(SC1 = clinical_data$SHEDDING_SC1[match(rowname, 
				  																					   table = clinical_data$CEL)],
				  			  SC2 = clinical_data$SYMPTOMATIC_SC2[match(rowname, 
				  																					   		  table = clinical_data$CEL)],
				  			  SC3 = clinical_data$LOGSYMPTSCORE_SC3[match(rowname, 
				  																					   		table = clinical_data$CEL)],
				  			  TIME = clinical_data$TIME[match(rowname, table = clinical_data$CEL)]) %>%
				 gather(Gene, expression, -rowname, -SC1, -SC2, -SC3, -TIME)				   				
outcomes <- c("SC1", "SC2", "SC3") 
logLS <- lapply(outcomes, function(OUTCOME) {
	if(OUTCOME == "SC3") {
	lDF <- matDF1 %>%
			   dplyr::select_(.dots = c("rowname", OUTCOME, "TIME", "Gene", "expression"))
	colnames(lDF)[colnames(lDF) %in% OUTCOME] = "outcome"
	lDF$expression[lDF$expression == 0] <- 0.25
	lDF <- lDF %>%
			   group_by(Gene, TIME) %>%
			   mutate(FC_sCor = cor.test(expression, outcome, method="spearman")$estimate) %>%
			   as.data.frame() %>%
			   dplyr::select(Gene, FC_sCor, TIME) %>%
			   unique() %>%
			   mutate(Outcome = OUTCOME)		
	} else {
	lDF <- matDF1 %>%
			   dplyr::select_(.dots = c("rowname", OUTCOME, "TIME", "Gene", "expression"))
	colnames(lDF)[colnames(lDF) %in% OUTCOME] = "outcome"
	lDF$expression[lDF$expression == 0] <- 0.25
	lDF <- lDF %>%
			   group_by(Gene, TIME) %>%
			   mutate(FC_sCor = log2(mean(expression[outcome %in% 1]) /
			   										  mean(expression[outcome %in% 0]))) %>%
			   as.data.frame() %>%
			   dplyr::select(Gene, FC_sCor, TIME) %>%
			   unique() %>%
			   mutate(Outcome = OUTCOME)
			}
		return(value = lDF)
		})
logDF <- do.call(what = rbind, logLS) %>%
				spread(Outcome, FC_sCor)


# load predictor lists files				   									   
fileLS <- list.files(pattern = "*Time.+_V3.xlsx$",
							path = ".",
							recursive = TRUE,
							full.names = TRUE)
							
# merge predictor lists of teams with their gene symbols
predLS <- lapply(fileLS, function(FILE) {
					print(FILE)
					datDF <- read.xls(FILE)
					datDF2 <- merge(annotDF, datDF, by="PROBEID", all.y = TRUE) %>%
									 mutate(GENESYMBOL = ifelse(GENE.y == "",
									 													GENE.x,
									 													GENE.y)) %>%
									 dplyr::select(-GENE.x, -GENE.y)
					return(value = datDF2)
				})
predDF <- do.call(rbind, predLS)


##### perform 'superExcatTest' to test significance of intersection of predictors between teams
challenges <- unique(predDF$CHALLENGE)
times <- unique(predDF$TIME)
SeTLS <- lapply(challenges, function(ch){
				lapply(times, function(tp) {
					print(ch)
					print(tp)
					teams <- read.xlsx("leaderboard.xlsx", sheetIndex = 1) %>%
			   					   filter(Challenge %in% ch & Time %in% tp ) %>%
			   					   filter(!ToUse %in% "No") %>%
			   					  .$team %>% as.character() %>% unique(.)	   
			   		# extract predictors for each timepoint and challenge
			   		genes <- logDF %>%
			   					   select_(.dots = c("Gene", "TIME", ch)) %>%
			   					   filter(TIME %in% tp) %>%
			   					   .$Gene 
					subPred <- predDF %>%
									   filter(CHALLENGE %in% ch & 
									   		   TIME %in% tp &
									   		   TEAM %in% teams & 
									   		   GENESYMBOL %in% genes)	%>%
									   dplyr::select(GENESYMBOL, TEAM)					   
				   teamLS <- unstack(subPred)
				   teamLS <- lapply(teamLS, unique)
	 		   # perform 'superExcatTest'
#	 		    seT <- MSET(x = teamLS,
	# 		   						  n = length(unique(annotDF$GENE)),
	 #		   					 	  lower.tail = FALSE)
	 			SeT <- supertest(teamLS, n = length(unique(annotDF$GENE)))
	 			summSeT <- summary(SeT)
	 			intersectionDF <- data.frame(pValue = summSeT$P.value,
	 														   ovSize = summSeT$otab) %>%
	 										  rownames_to_column() %>%
	 										  dplyr::rename(Sets = rowname) %>%
	 										  #filter(pValue < 0.05) %>%
	 										  mutate(Sets1 = strsplit(Sets, "")) %>%
	 										  group_by(Sets) %>%
	 										  mutate(nSets = sum(as.numeric(unlist(Sets1)))) %>%
	 										  as.data.frame() %>%
	 										  dplyr::select(-Sets1) %>%	 										  
	 										  group_by(nSets) %>%
	 										  mutate(nSigCombinations = as.numeric(table(pValue<0.05)["TRUE"]),
	 										  			  TotalCombsPerSet = as.numeric(length(Sets)),
	 										  			  perSigCombinations = (nSigCombinations/TotalCombsPerSet)*100) %>%
	 										  as.data.frame() %>%
	 										  group_by(nSets) %>%
	 										  mutate(medianOverlapSizeAmongSigCombs = 
	 										  			  median(ovSize[pValue<0.05])) %>%
	 										  as.data.frame() %>%
	 										  dplyr::select(nSets, 
	 										  					  nSigCombinations, 
	 										  					  TotalCombsPerSet,
	 										  					  perSigCombinations,
	 										  					  medianOverlapSizeAmongSigCombs) %>%
	 										  unique() %>%
	 										  arrange(nSets) %>%
	 										  mutate(CHALLENGE = ch,
	 										  			  TIMEPOINT = tp) %>%
	 										  dplyr::select(nSets, 
	 										  			 		  perSigCombinations, 
	 										  			 		  medianOverlapSizeAmongSigCombs,
	 										  			 		  CHALLENGE,
	 										  			 		  TIMEPOINT) %>%
	 										  filter(!nSets %in% 1)
	 		  return(value = intersectionDF)
	 		  })
	 	})
SeTDF <- do.call(rbind, lapply(SeTLS, function(x) do.call(rbind, x)))


### Super Test on SC2 TIME24
seLS <- list.files(path = "SuperTest_SC2_TIME24",
						   pattern = "*.RData",
						   recursive = FALSE,
						   full.names = TRUE)
seLS1 <- lapply(seLS, function(FILE) {
	rm(summSeT)
	load(FILE)
	intDF  <- data.frame(pValue = summSeT$P.value,
	 								 ovSize = summSeT$otab) 
	return(value = intDF)
			})
			
df2 <- do.call(rbind, seLS1) %>%
											 as.data.frame() %>%
	 										  rownames_to_column() %>%
	 										  dplyr::rename(Sets = rowname) %>%
	 										  #filter(pValue < 0.05) %>%
	 										  mutate(Sets1 = strsplit(Sets, "")) %>%
	 										  group_by(Sets) %>%
	 										  mutate(nSets = sum(as.numeric(unlist(Sets1)))) %>%
	 										  as.data.frame() %>%
	 										  dplyr::select(-Sets1) %>%	 										  
	 										  group_by(nSets) %>%
	 										  mutate(nSigCombinations = as.numeric(table(pValue<0.05)["TRUE"]),
	 										  			  TotalCombsPerSet = as.numeric(length(Sets)),
	 										  			  perSigCombinations = (nSigCombinations/TotalCombsPerSet)*100) %>%
	 										  as.data.frame() %>%
	 										  group_by(nSets) %>%
	 										  mutate(medianOverlapSizeAmongSigCombs = 
	 										  			  median(ovSize[pValue<0.05])) %>%
	 										  as.data.frame() %>%
	 										  dplyr::select(nSets, 
	 										  					  nSigCombinations, 
	 										  					  TotalCombsPerSet,
	 										  					  perSigCombinations,
	 										  					  medianOverlapSizeAmongSigCombs) %>%
	 										  unique() %>%
	 										  arrange(nSets) %>%
	 										  mutate(CHALLENGE = "SC2",
	 										  			  TIMEPOINT = "TIME24") %>%
	 										  dplyr::select(nSets, 
	 										  			 		  perSigCombinations, 
	 										  			 		  medianOverlapSizeAmongSigCombs,
	 										  			 		  CHALLENGE,
	 										  			 		  TIMEPOINT) %>%
	 										  filter(!nSets %in% 1)
plotDF <- rbind(df1, df2, df3, df4)


# bar plot
plotDF <- df4 %>% 
				mutate(medianOverlapSizeAmongSigCombs = 
							round(medianOverlapSizeAmongSigCombs))
outFile <- "PredictorIntersections_SC3_TIME24.pdf"
order <- plotDF$nSets
		plotTheme <- theme(panel.grid = element_blank(),
										 axis.text.x = element_text(size = 8, color = "black"),
										 axis.text.y = element_text(size = 8, color = "black"),
										 axis.title.x = element_text(size = 8, color = "black"),
										 axis.title.y = element_text(size = 8, color = "black"))
	  bP <- ggplot(data = plotDF,
	  					  mapping = aes(x = nSets, 
	  					  						   y = perSigCombinations, 
	  					  						   size = medianOverlapSizeAmongSigCombs)) +
	  			geom_point()  +
	  			scale_size_area(max_size = 10, 
	  									   breaks = sort(unique(plotDF$medianOverlapSizeAmongSigCombs)),
	  									   name = "Median intersection size") +
			   scale_x_discrete(limits = order, labels = order) +
			   labs(x = "Number of teams",
			   		   y = "% Team combinations showing a significant intersection of predictors") +
			   theme_bw() +
			   plotTheme
		pdf(file = file.path(resDir,outFile), width = 5, height = 3)
		print(bP)
		dev.off()
				





### Pathway enrichment (Fisher's exact Test)
# define geneSet background
  gmx <- "h.all.v6.0.symbols.gmt"
  dirName <- "DB"
  gmxFile <- file.path(dirName, gmx)
  colNames <- max(count.fields(file = gmxFile, sep = "\t"))
  colNames <- seq(from = 1, to = colNames)
  colNames <- as.character(colNames)
  gmx <- read.table(file      = gmxFile,
                    sep       = "\t",
                    quote     = "\"",
                    fill      = TRUE,
                    col.names = colNames,
                    row.names = 1)
 gmx <- gmx[, -1]
 gmx <- apply(gmx, MARGIN = 1, FUN = function(x) {
    						return(value = setdiff(unname(x), ""))
  					})
 names(gmx) <- toupper(names(gmx))

 # extract the name of db from the gmx arguments
 db <- basename(gmxFile)
 db <- gsub(pattern = "[.].+$",
            		replacement = "",
            		db)
 # obtain background
 bg <- unique(unlist(gmx))
 fisherIndex <- NULL


##### pathway enrichment only across top scoring teams (leaderboard > 50% AUC)
challenges <- unique(predDF$CHALLENGE)
times <- unique(predDF$TIME)
resDir = "results_20171102"
fisherLS <- lapply(challenges, function(ch) {
						lapply(times, function(tp) {
					print(ch)
					print(tp)
					teams <- read.xlsx("leaderboard.xlsx", sheetIndex = 1) %>%
			   					   filter(Challenge %in% ch & Time %in% tp ) %>%
			   					   filter(!ToUse %in% "No") %>%
			   					  .$team %>% as.character() %>% unique(.)	   
			   		# extract predictors (positive FC / sCor for the outcome)	
			   		genes <- logDF %>%
			   					   select_(.dots = c("Gene", "TIME", ch)) %>%
			   					   filter(TIME %in% tp)
			   		colnames(genes)[colnames(genes) %in% ch] <- "outcome"
			   		genesP <- genes %>% 
			   					     filter(outcome > 0) %>%
			   					     .$Gene 
			   		genesN <- genes %>% 
			   					     filter(outcome < 0) %>%
			   					     .$Gene 	
			   		# subset team predictors on the genes above	   
					subPredP <- predDF %>%
									     filter(CHALLENGE %in% ch & 
									   		   TIME %in% tp &
									   		   TEAM %in% teams & 
									   		   GENESYMBOL %in% genesP) %>%
									   	 mutate(ASSOCIATION = "Pathways associated with symptoms")						   
					subPredN <- predDF %>%
									     filter(CHALLENGE %in% ch & 
									   		   TIME %in% tp &
									   		   TEAM %in% teams & 
									   		   GENESYMBOL %in% genesN) %>%
									   	 mutate(ASSOCIATION = "Pathways associated with lack of symptoms")					   
				   subPred <- rbind(subPredP, subPredN)
				   subPredLS <- lapply(unique(subPred$ASSOCIATION), function(a) {
					gs <- subPred %>%
							 filter(ASSOCIATION %in% a) %>%
		 					 .$GENESYMBOL %>%
		 					  unique(.)
							output <- mclapply(gmx, function(x) {
											 tab <- table(factor(bg %in% gs, levels = c(TRUE, FALSE)),
             													factor(bg %in% x, levels = c(TRUE, FALSE)))
							fit <- fisher.test(tab, alternative = "greater")
							interSection <- intersect(gs, x)
							interSection <- paste(interSection, collapse = ",")
							return(value = c(RATIO = as.numeric((diag(tab) / colSums(tab))[1]),
                 									  NOM_pval = fit$p.value,
                 									  INTERSECT = interSection))
       												})
						    fDF <- do.call(what = rbind, args = output) %>%
						    		   		as.data.frame() %>%
						    		   		mutate(NAME = names(gmx),
						    			   		   		ADJ_pval = p.adjust(as.numeric(NOM_pval), 
						    			   			   									 method = "BH"),
						    			   			   	ASSOCIATION = a) 
							})
					subPredDF <- do.call(rbind, subPredLS)
					# save results
					 outFile <- subPredDF %>% 
					 				 filter(ADJ_pval < 0.05) %>% 
					 				 arrange(ADJ_pval) %>%
						    		 dplyr::select(NAME, RATIO, NOM_pval, ADJ_pval, INTERSECT, ASSOCIATION) %>%
						    		     mutate(CHALLENGE = rep(ch, length(NAME)),
						    		   			     TIME = rep(tp, length(NAME)),
						    		   	 		     log10P = -log10(ADJ_pval))
					 outFileName <- paste(ch, tp, "PathwaysAcrossTeams_TOPLB_PosandNeg.txt", sep="_")
					 write.table(outFile, 
					 				   file = file.path(resDir, outFileName),
					 				   quote = FALSE,
					 				   sep = "\t",
					 				   row.names = FALSE) 
					return(value = outFile)
				})
			})
fisherDF <- do.call(rbind, lapply(fisherLS, function(x) do.call(rbind, x)))

# plot pathways
		plotDF <- fisherDF %>%
						mutate(NAME = gsub("HALLMARK_", "", NAME)) %>%
						group_by(NAME, CHALLENGE, TIME) %>%
						mutate(mostSig = max(log10P)) %>%
						as.data.frame() %>%
						mutate(mostSigAssociation = ifelse(log10P==mostSig, "keep", "no")) %>%
						filter(mostSigAssociation %in% "keep")
# plot
		plotTheme <- theme(panel.grid = element_blank(),
										 axis.text.x = element_text(size = 8, color = "black"),
										 axis.text.y = element_text(size = 8, color = "black"),
										 axis.title.x = element_text(size = 8, color = "black"),
										 axis.title.y = element_text(size = 8, color = "black"))
		cols <- c("Pathways associated with symptoms" = "red", 
					  "Pathways associated with lack of symptoms" = "blue")

# order by pathway with maximum occurance
	orderP <- plotDF %>% dplyr::select(NAME, ASSOCIATION) %>%
				    arrange(ASSOCIATION) %>%
				    unique() %>%
				    filter(ASSOCIATION %in% "Pathways associated with symptoms") %>%
				    .$NAME %>%
				    rev()
	orderN <- plotDF %>% dplyr::select(NAME, ASSOCIATION) %>%
				    arrange(ASSOCIATION) %>%
				    unique() %>%
				    filter(ASSOCIATION %in% "Pathways associated with lack of symptoms") %>%
				    .$NAME %>%
				    rev()
order <- c(orderN, orderP)
# bubble plot
	figDir = "results_20171102"
	outFile <- "Pathways across top scoring teams on Leaderboard_PositiveandNegative_V2.pdf"
	  bP <- ggplot(plotDF, aes(CHALLENGE, NAME, color = ASSOCIATION, size = log10P)) + 
				geom_point() +
				scale_size(range = c(1,10), name = "-log10(p-value)") +
			    scale_color_manual(values = cols, name = NULL) +
			   facet_wrap(~TIME) +
			   scale_y_discrete(limits = order) +
			   labs(y="PATHWAY") +
			   theme_bw() +
			   plotTheme
		pdf(file = file.path(figDir,outFile), width = 8, height = 7)
		print(bP)
		dev.off()



# list results
resLS <- list.files(pattern = "*_TOPLB_Positive.txt",
							path = resDir,
							full.names = TRUE,
							recursive = FALSE)

# features annotation
featDF <- data.frame(cName = colnames(expression_data_gene_level)) %>%
				mutate(probeID = gsub(pattern = "(.+):.+", replacement = "\\1", cName),
							geneSymbol = gsub(pattern = ".+:(.+)", replacement = "\\1", cName))

# expression data
eDat <- expression_data_gene_level
colnames(eDat) <- gsub(pattern = "(.+):.+", replacement = "\\1", colnames(eDat))


# do SLEA (Sample Level Enrichment Analysis) function(z-score per sample)
doSLEA <- function(expressionSet, geneSet) {
			# scale expression
						exprsMat  <- expressionSet
						exprsMat  <- t(scale(t(exprsMat)))
			# extract expression of leGenes of each geneset
						comm <- intersect(geneSet, rownames(expressionSet))
	     				gsDF <- exprsMat[comm, ]
			# calculate mean expression per sample
						gsM <- colMeans(gsDF)	
			# extract random genes of size of the geneSet from full probeset and calculate mean
			# and perform this for 'n' permutations
						nperm <- lapply(1:1000, function(j) {
			# set seed for every permutation
	    												set.seed(j)
	    												rGSDF <- exprsMat[sample.int(nrow(exprsMat),
	    																									length(comm)), ] 
	    												rGSM <- colMeans(rGSDF)
													return(value = rGSM)
												})
						permDF <- do.call(rbind, nperm)	
						zscore <- (gsM - colMeans(permDF)) / apply(permDF,2,sd) 	
						sleaDF <- zscore %>% as.data.frame()	  
						return(value = sleaDF)
					}


# Iterate through every result
for(RES in resLS) {
	print(RES)
	resDF <- read.delim(file = RES, stringsAsFactors = FALSE)
	challenge <- gsub(pattern = ".+(SC.+)_TIME.+", replacement = "\\1", RES)
	time <- gsub(pattern = ".+SC.+_(TIME.+)_Pathways.+", replacement = "\\1", RES)
		#for(i in c(1:nrow(resDF))) {
			lapply(1:nrow(resDF), function(i) {
					if(grepl("SC1|SC2", challenge)) {
			pathway <- resDF[i, ] %>% .$NAME %>% as.character()
			#print(pathway)
			genes <-  resDF %>% 
							filter(NAME %in% pathway) %>% 
							.$INTERSECT %>%
							strsplit(",") %>%
							unlist()
			probes <- featDF %>% filter(geneSymbol %in% genes) %>% .$probeID
	
	# subset samples on TIME0 / TIME24
	subSamp <- clinical_data
	outcome <- colnames(subSamp)[grep(challenge, colnames(subSamp))]
	subSamp <- subSamp %>% filter(TIME %in% time ) %>%
						 dplyr::select_(.dots = c("STUDYID", 
						 							#"SUBJECTID", 
						 							"AGE", 
						 							"GENDER",
						 							outcome,
						 							"TIMEHOURS",
						 							"CEL",
						 							"VIRUS"))
	rownames(subSamp) <- subSamp$CEL		
	outCol <- grep(challenge, colnames(subSamp))
	subSamp <- subSamp[order(subSamp[, outCol]), ]
	rownames(subSamp) <- subSamp$CEL
	colnames(subSamp)[outCol] <- challenge
	subSamp <- subSamp %>% dplyr::select(-CEL)

# call doSLEA
expressionSet <- eDat[rownames(subSamp), ] %>% t() %>% as.data.frame()
geneSet <- probes
sleaDF <- doSLEA(expressionSet = expressionSet, geneSet = geneSet)
names(sleaDF) = pathway
sleaDF <- sleaDF %>% 
				 rownames_to_column() %>%
				 dplyr::rename(CEL = rowname)
subSamp1 <- subSamp %>% rownames_to_column() %>%
					   dplyr::rename(CEL = rowname)
sleaDF <-  merge(sleaDF, subSamp1, by = "CEL", all.x = TRUE) %>% as.data.frame()

# define theme of bar plot
plotTheme <- theme(axis.title.y = element_text(size=12),
								 axis.ticks.x = element_blank(),
								 axis.text.x = element_blank(),
				   				 axis.text.y = element_text(size=12, color="black"),
				   				 panel.background = element_blank(),
				   				 panel.grid = element_blank(),
				   				 axis.line.x = element_line(color="black"),
				   				 axis.line.y = element_line(color="black"))
cols <- c("0" = "blue", "1"  ="red")
#order <- c(0, 1)
yName <- gsub("HALLMARK_", "", pathway)

# jitter plot of slea z-score
plotJit <- ggplot(data = sleaDF,
						  mapping = aes(x = as.factor(sleaDF[, challenge]), 
						  						   y = sleaDF[, pathway])) +			   
			   geom_boxplot(outlier.colour = NA, size = 2) +		 
			   geom_jitter(mapping = aes(color = as.factor(sleaDF[, challenge])), height = 0, size = 2) +			   
			   scale_color_manual(values = cols, name = challenge) +
			   #scale_x_discrete(limits = order) +
			   labs(x = NULL, y = yName) +
				plotTheme
 outFile2 <- paste(pathway,
	  						 challenge, time, "SLEAzscore_TOPLB_Positive.pdf", sep= "_")   		    
 pdf(file = file.path(figDir, outFile2), width = 5, height = 4)
 print(plotJit)
 dev.off()

# For SC1 and SC2 : wilcox tests
print(paste(pathway, challenge, time, sep=" : "))
wT <- wilcox.test(sleaDF[, pathway] ~ as.factor(sleaDF[, challenge]))
print(wT$p.value)
 } else {
# for SC3
		pathway <- resDF[i, ] %>% .$NAME %>% as.character()
		genes <-  resDF %>% 
							filter(NAME %in% pathway) %>% 
							.$INTERSECT %>%
							strsplit(",") %>%
							unlist()
			probes <- featDF %>% filter(geneSymbol %in% genes) %>% .$probeID
	
	# subset samples on TIME0 / TIME24
	subSamp <- clinical_data
	outcome <- colnames(subSamp)[grep(challenge, colnames(subSamp))]
	subSamp <- subSamp %>% filter(TIME %in% time ) %>%
						 dplyr::select_(.dots = c("STUDYID", 
						 							#"SUBJECTID", 
						 							"AGE", 
						 							"GENDER",
						 							outcome,
						 							"TIMEHOURS",
						 							"CEL",
						 							"VIRUS"))
	rownames(subSamp) <- subSamp$CEL		
	outCol <- grep(challenge, colnames(subSamp))
	subSamp <- subSamp[order(subSamp[, outCol]), ]
	rownames(subSamp) <- subSamp$CEL
	colnames(subSamp)[outCol] <- challenge
	subSamp <- subSamp %>% arrange(SC3)
	rownames(subSamp) <- subSamp$CEL
	subSamp <- subSamp %>% dplyr::select(-CEL)

# call doSLEA
		expressionSet <- eDat[rownames(subSamp), ] %>% t() %>% as.data.frame()
		geneSet <- probes
		sleaDF <- doSLEA(expressionSet = expressionSet, geneSet = geneSet)
		names(sleaDF) = pathway
		sleaDF <- sleaDF %>% 
				 		rownames_to_column() %>%
				 		dplyr::rename(CEL = rowname)
		subSamp1 <- subSamp %>% rownames_to_column() %>%
					   			dplyr::rename(CEL = rowname)
		sleaDF <-  merge(sleaDF, subSamp1, by = "CEL", all.x = TRUE) %>% as.data.frame()
		sleaDF2 <- sleaDF %>% arrange_(.dots = pathway)
		order <- sleaDF2$CEL
		yName <- gsub("HALLMARK_", "", pathway)

# pathway correlation with outcome		
		print(paste(pathway, challenge, time, sep=" : "))
		cT <- cor.test(sleaDF2[, 2], sleaDF2[, challenge], method="spearman")
		print(paste(round(cT$p.value, 4), round(cT$estimate, 3), sep = " : "))
	# bar plot of pathway z-score
		plotbar <- ggplot(data  = sleaDF2,
						   	mapping = aes(x = CEL, y = sleaDF2[, 2], fill = SC3)) +
						   	scale_fill_gradient(low = "blue", high = "red") +
						   	geom_bar(stat = "identity") +
				  		    labs(x = NULL, y = yName) +
				   	    		scale_x_discrete(limits = order) +
				   		    plotTheme
 		outFile2 <- paste(pathway,
	  								 challenge, time, "SLEAzscore_TOPLB_Positive.pdf", sep= "_")    
 		pdf(file = file.path(figDir, outFile2), width = 15, height = 4)
 		print(plotbar)
 		dev.off()
	}
})
}
 
 
# check for significance of gene set overlap between SC2 and SC3 at Time0 and Time24
res2LS <- list.files(path = resDir,
							  pattern = "_Positive.txt$",
							  full.names = TRUE,
							  recursive = FALSE)
res2LS <- res2LS[grep("SC2|SC3", res2LS)]
tps <- c("TIME0", "TIME24")
for(TP in tps) {
#ovLS <- lapply(tps, function(i) {
	print(TP)
			fileLS <- res2LS[grep(TP, res2LS)]
		# pathways common between SC2 and SC3
			comm <- intersect(read.delim(fileLS[1]) %>% .$NAME,
										  read.delim(fileLS[2]) %>% .$NAME)
			
	# check overlaps							
for(PATHWAY in comm) {
	sc2 <- read.delim(fileLS[1]) %>%
			   filter(NAME %in% PATHWAY) %>%
			   .$INTERSECT %>%
				strsplit(",") %>%
				unlist(.)
	sc3 <- read.delim(fileLS[2]) %>%
			   filter(NAME %in% PATHWAY) %>%
			   .$INTERSECT %>%
				strsplit(",") %>%
				unlist(.)
# make contingency table
	conTab <- table(factor(gmx[[PATHWAY]] %in% sc2, levels = c(TRUE,FALSE)),
                			  factor(gmx[[PATHWAY]] %in% sc3, levels = c(TRUE,FALSE)))
    fisherT <- fisher.test(conTab, alternative = "greater")
	expOvl <- (length(sc2) * length(sc3)) / length(gmx[[PATHWAY]])
	obsOvl <- conTab[1,1]		
	if(fisherT$p.value <= 0.05) {
		  if(obsOvl > expOvl) {
				message <- paste(PATHWAY, obsOvl, expOvl, fisherT$p.value,
											  "obsOvl > expOvl : There is a significant overlap of genes",
											   sep = " : ")
				print(message)
					} else {
				message <- paste(PATHWAY, obsOvl, expOvl,fisherT$p.value,
											  "obsOvl < expOvl : There is a significant difference in the leGenes", sep = " : ")
						print(message)
					}
				} else {
					message <- paste(PATHWAY,
												" : Fisher test is not significant", sep = "")
					print(message)
					print(obsOvl)
					print(expOvl)
				}
		}
	}




# download genemani JAR file
genemaniaPath <- file.path(".", "plugin-cy3-3.4.1.jar")
genemaniaFile <- basename(genemaniaPath)

# identify most recent genemania database
genemaniaJava <- "java -Xmx1G -jar"
genemaniaCmd <- "DataAdmin list"
genemaniaCall <- paste(genemaniaJava, genemaniaFile, genemaniaCmd)
genemaniaIntern <- system(command  = genemaniaCall,
                          					 intern  = TRUE,
                          					 ignore.stderr = TRUE)
genemaniaOutput <- genemaniaIntern %>%
  								 strsplit(split = "\t") %>%
  								 do.call(what = rbind)
colnames(genemaniaOutput) <- genemaniaOutput[1, ]
genemaniaDB <- genemaniaOutput %>%
  							as.data.frame() %>%
  							filter(grepl(pattern = "core", `Data Set ID`)) %>%
  							arrange(desc(`Data Set ID`)) %>%
  							.[1, "Data Set ID"]

# download genemania database
genemaniaCmd <- "DataAdmin install"
genemaniaCall <- paste(genemaniaJava, genemaniaFile, genemaniaCmd, genemaniaDB)
genemaniaIntern <- system(command  = genemaniaCall,
                         					 intern  = TRUE,
                         					 ignore.stderr = TRUE)
genemaniaDB <- paste0("gmdata-", genemaniaDB)
genemaniaCmd <- "DataAdmin list-data"
genemaniaCall <- paste(genemaniaJava, genemaniaFile, genemaniaCmd, genemaniaDB)
genemaniaIntern <- system(command  = genemaniaCall,
                         	 				 intern  = TRUE,
                          					 ignore.stderr = TRUE)
genemaniaOutput <- genemaniaIntern %>%
  								  gsub(pattern = "\t$", replacement = "\t ") %>%
  								  strsplit(split = "\t") %>%
  								  do.call(what = rbind)
colnames(genemaniaOutput) <- genemaniaOutput[1, ]
genemaniaID <- genemaniaOutput %>%
  						   as.data.frame() %>%
  						   filter(grepl(pattern = "Human", Description)) %>%
  						   .$"Data ID"
genemaniaCmd <- "DataAdmin install-data"
genemaniaCall <- paste(genemaniaJava,
                       				   genemaniaFile,
                       				   genemaniaCmd,
                       				   genemaniaDB,
                       				   genemaniaID)
genemaniaIntern <- system(command = genemaniaCall,
                          					intern  = TRUE,
                          					ignore.stderr = TRUE)


# GeneMania based on leading edge genes
# prepare input file
genemaniaParamsPath <- file.path(".", "genemania-parameters.csv")
genemaniaDefault <- scan(file = genemaniaParamsPath, what = "raw", sep = "\n")
genemaniaJava <- "java -Xmx1G -jar"
genemaniaJar <- list.files(pattern = "plugin.+jar")
genemaniaCmd <- "QueryRunner"
genemaniaDB <- list.files(pattern = "gmdata", include.dirs = TRUE)
genemaniaDir <- "geneMania"
queryFile <- "genemania_genesig.query"
logFileName <- gsub(pattern = ".query", replacement = ".log", queryFile)
relatedGeneLimit <- 0


# write input file
# Gene mania network of genes
ovLS <- lapply(tps, function(TP) {
	print(TP)
			fileLS <- res2LS[grep(TP, res2LS)]
   # pathways common between SC2 and SC3
			comm <- intersect(read.delim(fileLS[1]) %>% .$NAME,
										  read.delim(fileLS[2]) %>% .$NAME)
			
	# check overlaps							
	lapply(comm, function(PATHWAY) {
	print(PATHWAY)
	sc2 <- read.delim(fileLS[1]) %>%
			   filter(NAME %in% PATHWAY) %>%
			   .$INTERSECT %>%
				strsplit(",") %>%
				unlist(.)
	sc3 <- read.delim(fileLS[2]) %>%
			   filter(NAME %in% PATHWAY) %>%
			   .$INTERSECT %>%
				strsplit(",") %>%
				unlist(.)
	commGenes <- intersect(sc2, sc3)
	pDF <- data.frame(commGenes = commGenes,
								  Pathway = PATHWAY,
								  Time = TP)
	})
})
ovDF <- do.call(rbind, lapply(ovLS, function(x) do.call(rbind, x)))


#### Table of number of teams that select a predictor

# remove the teams not to use from predDF table
keepTeams <- read.xlsx("leaderboard.xlsx", sheetIndex = 1) %>%
					   mutate(TChT = interaction(team, Challenge, Time, sep = "_")) %>%
					   filter(!ToUse %in% "No") %>%
					  .$TChT %>%
					   unique(.) %>%
					   as.character()
teamDF <- predDF %>%
				  mutate(TChT = interaction(TEAM, CHALLENGE, TIME, sep = "_")) %>%
				  filter(TChT %in% keepTeams & !CHALLENGE %in% "SC1") %>%
				  filter(GENESYMBOL %in% unique(ovDF$commGenes)) %>%
				  group_by(GENESYMBOL, TIME) %>%
				  mutate(nTeams = length(unique(TEAM))) %>%
				  as.data.frame() %>%
				  dplyr::select(TIME, GENESYMBOL, nTeams) %>%
				  unique()


## Plot network of TIME0 and TIME24 genes
timePoint <- "TIME24"

if(timePoint == "TIME0") {
gs <- ovDF %>%
	     filter(Time %in% timePoint & Pathway %in% c("HALLMARK_HEME_METABOLISM", 
	     																			"HALLMARK_INFLAMMATORY_RESPONSE",
	     																			"HALLMARK_KRAS_SIGNALING_UP")) %>%
	     .$commGenes %>%
	     unique(.)
} else {
	gs <- ovDF %>%
			filter(Time %in% timePoint) %>%
		     filter(Pathway %in% c("HALLMARK_HEME_METABOLISM",
						  						 "HALLMARK_INFLAMMATORY_RESPONSE",
						  						 "HALLMARK_KRAS_SIGNALING_UP",
						  						 "HALLMARK_INTERFERON_GAMMA_RESPONSE")) %>%
	     	.$commGenes %>%
	     	unique(.)
	}	     
	
  write("H. Sapiens", file = queryFile)
  write(paste(gs, collapse = "\t"), file = queryFile, append = TRUE)
  write(paste(c("coexp", "pi"), collapse = "\t"), file = queryFile, append = TRUE)  
  write(paste(genemaniaDefault, collapse = "\t"),
        			file   = queryFile,
        			append = TRUE)
  write(relatedGeneLimit, file = queryFile, append = TRUE)
  write("average", file = queryFile, append = TRUE)

 # call genemania
  genemaniaCall <- paste(genemaniaJava,
                         				 genemaniaJar,
                         				 genemaniaCmd,
                         				queryFile,
                         				"--data", genemaniaDB,
                         				"--out xml",
                         				"--results",
                         				genemaniaDir,
                         				">",
                         				logFileName)
  genemaniaIntern <- system(command = genemaniaCall,
                            				   intern = TRUE,
                            				   ignore.stderr = TRUE)

 outputFile <- paste0(queryFile, "-results.report.xml")
 outputFile <- file.path(genemaniaDir, basename(outputFile))
 textXML <- xmlRoot(xmlTreeParse(file = outputFile))
 textDF <- xmlChildren(textXML[["results"]][["interactions"]]) %>%
   				 lapply(FUN = xmlAttrs) %>%
    			 do.call(what = rbind) %>%
    			 as.data.frame(row.names = as.character(1:nrow(.))) %>%
    			 mutate(weight = as.numeric(weight))
 symbol <- union(textDF$from, textDF$to)

# define edge attributes
  edgesDF <- textDF %>%
    				  dplyr::select(from, to, type) %>%
    				  distinct() %>%
    				  mutate(color = c("Co-expression" = "thistle",
                     							"Physical Interactions" = "green")[type])
# make graph
  g <- graph.empty(directed = FALSE) %>%
    	 add.vertices(nv = length(symbol), name = symbol, color = NA) %>%
      	 add.edges(edges  = t(edgesDF[, c("from", "to")]),
                		   weight = edgesDF$weight,
                		   type   = edgesDF$type,
                		   color  = edgesDF$color)

# define node sizes
if(timePoint == "TIME0") {
			nodeDF <- teamDF %>% 
							  filter(TIME %in% timePoint & GENESYMBOL %in% V(g)$name)
			nodeDF <- nodeDF[order(match(nodeDF$GENESYMBOL, 
																table = V(g)$name)), ]	
									} else {
			nodeDF <- teamDF %>%
							  filter(TIME %in% timePoint & GENESYMBOL %in% V(g)$name)
			nodeDF <- nodeDF[order(match(nodeDF$GENESYMBOL, 
																 table = V(g)$name)), ]
								  }
maxSize <- max(nodeDF$nTeams)

# define vertex attributes for tk plot positioning
if(timePoint == "TIME0") {
	sub_ovDF <- ovDF %>% filter(Time %in% timePoint) %>%
		    			  filter(Pathway %in% c("HALLMARK_HEME_METABOLISM",
						  						 			  "HALLMARK_INFLAMMATORY_RESPONSE",
						  						 			  "HALLMARK_KRAS_SIGNALING_UP")) %>%
						  mutate(color = ifelse(Pathway %in% "HALLMARK_HEME_METABOLISM",
						  									"deepskyblue2",
						  						  ifelse(Pathway %in% "HALLMARK_INFLAMMATORY_RESPONSE",
						  						  			"orange", "indianred1")))
		} else {
	sub_ovDF <- ovDF %>% filter(Time %in% timePoint) %>%
		    			  filter(Pathway %in% c("HALLMARK_HEME_METABOLISM",
						  						 			  "HALLMARK_INFLAMMATORY_RESPONSE",
						  						 			  "HALLMARK_KRAS_SIGNALING_UP",
						  						 			  "HALLMARK_INTERFERON_GAMMA_RESPONSE")) %>%
						  mutate(color = ifelse(Pathway %in% "HALLMARK_HEME_METABOLISM",
						  									"deepskyblue2",
						  						  ifelse(Pathway %in% "HALLMARK_INFLAMMATORY_RESPONSE",
						  						  			"orange", 
						  						  ifelse(Pathway %in% "HALLMARK_KRAS_SIGNALING_UP", "indianred1", "red"))))
		}
# define vertex colors
	vertexDF <- data.frame(vName = V(g)$name) %>%
						mutate(color = sub_ovDF$color[match(vName, 
																					   table = sub_ovDF$commGenes)])
						
# tkplot to get coordinates
tkplot(g, vertex.color = vertexDF$color, vertex.size = nodeDF$nTeams/maxSize * 30)
posMat <- tkplot.getcoords(tkp.id = 3)

outFile <- paste(timePoint, "_GeneMania_Network.pdf", sep = "")
pdf(file = file.path(resDir, outFile),width = 7, height = 6)
plot(g, 
	     vertex.shape = "circle", 
	     vertex.size = nodeDF$nTeams/maxSize * 30, 
	     vertex.color = vertexDF$color,
	     edge.color = edgesDF$color,
	     vertex.frame.color = "darkgrey",
	     vertex.label.cex = 0.5, 
	      layout = posMat,
	     vertex.label.color = "black")
legend(x = -1,
            y = -1.1,
            pch = 16,
            col = unique(sub_ovDF$color),
            legend = unique(sub_ovDF$Pathway),
            cex = 0.5,
            bty = "n",
            xjust = 0.5,
            title = "Pathway")
legend(x = 0,
            y = -1.1,
            lwd = 2,
            col = unique(edgesDF$color),
            legend = unique(edgesDF$type),
            cex = 0.5,
            bty = "n",
            xjust = 0.5,
            title = "Edge Annotation")
legend(x = 1,
            y = -1.1,
            pch = 16,
            col = c("black", "black"),
            legend = signif(round(range(nodeDF$"nTeams"), 2), digits = 3),
            cex = 0.5,
            pt.cex = c(1, 2),
            bty = "n",
            xjust = 0.5,
            title = "Number of teams that selected the predictor")
dev.off()    

	





