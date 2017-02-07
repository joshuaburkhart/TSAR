print("TSAR: Predict 0h")
print(" author: Joshua Burkhart")
print(" date: August 28, 2016")

print("libs")
library(knitr)
library(magrittr)
library(dplyr)
library(affy)
library(ggplot2)
library(limma)
library(pheatmap)
library(snm)
library(e1071)

print("globals")
IMAGES_AS_PNG <- TRUE
DATA_DIR <- "~/Downloads/" #"/media/burkhart/D052-7853/McWeeney/TSAR/"
PHAS1_CLINICL <- paste(DATA_DIR,"RespiratoryViralChallenge_IndependentTest_Biochronicity_Time0_CLINICAL.tsv",sep="")
                       #"ViralChallenge_test_Phase1_CLINICAL.tsv",sep="")
PHAS1_EXPRESS <- paste(DATA_DIR,"RespiratoryViralChallenge_IndependentTest_Biochronicity_Time0_EXPRESSION_VST.tsv",sep="")
                       #"ViralChallenge_test_Phase1_EXPRESSION_RMA.tsv",sep="")

print("misc_config")
rnaSeq2GeneSymbol_df <- read.table("/Users/joshuaburkhart/SoftwareProjects/TSAR/misc/anno_refseq.csv",header = T, sep=",") %>%
  dplyr::mutate(gene_name = as.character(gene_name)) %>%
  dplyr::mutate(refseq = as.character(refseq)) %>%
  dplyr::select(gene_name,refseq)
geneSymbol2probesetId_df <- read.table("/Users/joshuaburkhart/SoftwareProjects/TSAR/misc/HG-U133A_2.na36.annot.csv",header = T,sep=",") %>%
  dplyr::mutate(gene_name = as.character(Gene.Symbol)) %>%
  dplyr::mutate(probesetId = as.character(Probe.Set.ID)) %>%
  dplyr::select(gene_name,probesetId)

rnaSeq2probesetId_df <- dplyr::inner_join(rnaSeq2GeneSymbol_df,geneSymbol2probesetId_df,by="gene_name") %>%
  dplyr::select(refseq,probesetId)


setwd(DATA_DIR)

print("load features")
load("/Users/joshuaburkhart/SoftwareProjects/TSAR/src/top_2000.rda")
  #"/home/burkhart/Software/TSAR/src/top_2000.rda") #loads top_2000 from train_0h

print("load trained model")
load("/Users/joshuaburkhart/SoftwareProjects/TSAR/src/tunedModel.rda")
  #"/home/burkhart/Software/TSAR/src/tunedModel.rda") #loads trained_svm from train_0h

print("load phas1 clinical data")
phas1_clinicl_df <- read.table(PHAS1_CLINICL, header = T, sep = "\t")

print("load phas1 eset")
exprs_df <- read.table(
  PHAS1_EXPRESS,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE,
  as.is = TRUE
)

refseq <- rownames(exprs_df)
rownames(exprs_df) <- NULL
exprs_df <- cbind(refseq,exprs_df)

exprs_df <- dplyr::inner_join(exprs_df,rnaSeq2probesetId_df,by="refseq")
exprs_df <- exprs_df %>%
  subset(!duplicated(probesetId))

median_val <- exprs_df %>%
  dplyr::select(-refseq,-probesetId) %>%
  as.matrix() %>%
  median()

featureProbes <- data.frame(probesetId = top_2000$probesetid)
exprs_df <- dplyr::left_join(featureProbes,exprs_df,by="probesetId")

rownames(exprs_df) <- exprs_df$probesetId
exprs_df <- exprs_df %>%
  dplyr::select(-refseq,-probesetId)

rma_exprs <- as.matrix(
  exprs_df
)

rma_exprs[is.na(rma_exprs)] <- median_val

print("keep only t<0 arrays")
tlt0_cel <- phas1_clinicl_df %>%
  dplyr::filter(TIMEHOURS <= 0) %>%
  dplyr::mutate(SAMPLEID = as.character(SAMPLEID)) %>%
  dplyr::select(SAMPLEID)

tlt0_cel <- tlt0_cel %>% filter(SAMPLEID %in% colnames(rma_exprs))
phas1_clinicl_df <- phas1_clinicl_df %>% filter(SAMPLEID %in% tlt0_cel$SAMPLEID)

tlt0_exprs <- rma_exprs[, tlt0_cel$SAMPLEID]

rma_eset <- ExpressionSet(assayData = tlt0_exprs)

print("keep only training features")
tlt0_exprs <- rma_eset[top_2000$probesetid, ] %>%
  exprs() %>%
  as.data.frame() %>%
  dplyr::select(match(tlt0_cel %>%
                        dplyr::select(SAMPLEID) %>% .[,1],names(.)))

phas1_data <- tlt0_exprs %>%
  t() %>%
  as.data.frame() %>%
  dplyr::add_rownames(var = "SAMPLEID") %>%
  dplyr::select(-SAMPLEID)

print("predict phas1 scores")
phas1_predictions <- predict(tunedModel,phas1_data) %>% as.numeric()

leaderboard_submission_0h <- data.frame(
  LOGSYMPTSCORE_SC3 = log10(phas1_predictions)) %>%
  cbind(phas1_clinicl_df %>%
          dplyr::filter(TIMEHOURS <= 0) %>%
          dplyr::select(SUBJECTID),.)

print("store phas1 predictions")
write.table(leaderboard_submission_0h,
            file="/Users/joshuaburkhart/SoftwareProjects/TSAR/src/burkhart_0h_IndependentTest_Biochronicity_SC3_predictions.csv",
            #"/home/burkhart/Software/TSAR/src/burkhart_0h_SC3_predictions.csv",
            sep=",",
            row.names=FALSE,
            quote=FALSE)


