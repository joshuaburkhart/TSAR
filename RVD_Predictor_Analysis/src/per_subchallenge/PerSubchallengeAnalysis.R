# ---
# title: Per Timepoint Analysis
# author: Joshua Burkhart
# date: Aug 20, 2017
# ---

# load libraries
library(SuperExactTest)
library(org.Hs.eg.db)
library(VennDiagram)
library(hgu133a2.db)
library(ReactomePA)
library(annotate)
library(magrittr)
library(ggplot2)
library(dplyr)

# helper functions
ReactomeBarplot <- function(query_gene_symbols,bar_color,plot_title){
  
  entrez_query <- AnnotationDbi::select(org.Hs.eg.db, query_gene_symbols, "ENTREZID", "SYMBOL")
  
  entrez_query <- entrez_query[!duplicated(entrez_query[,1]), 2]
  
  reactomeEnrich <- ReactomePA::enrichPathway(gene=unique(entrez_query),
                                              organism = "human",
                                              pvalueCutoff = 1.0,
                                              pAdjustMethod = "BH",
                                              qvalueCutoff = 1.0,
                                              universe = ENTREZ_UNIVERSE,
                                              minGSSize = 2,
                                              maxGSSize = 1000,
                                              readable = T)
  reactomeEnrichAtLeast2GenesPerPathway <- reactomeEnrich %>% 
    as.data.frame()
  reactomeEnrichSummary <- reactomeEnrich %>% summary()
  
  qc_df <- data.frame(
    Pathway=as.character(reactomeEnrichAtLeast2GenesPerPathway$Description),
    log10.p.value=log10(reactomeEnrichAtLeast2GenesPerPathway$pvalue)
  )
  
  qc_df %>% dplyr::slice(1:30) %>% ggplot(
    aes(x=reorder(Pathway,-log10.p.value),y=-log10.p.value)) +
    theme(legend.position="none") +
    geom_bar(stat="identity", fill=bar_color) +
    geom_text(position="stack",aes(label=round(-log10.p.value,digits=3),hjust=1.1)) +
    coord_flip() +
    labs(x="Reactome Pathway",y="-log10(p-value)",title=plot_title)
}

## "an element is omitted if it is equal to any previous element"
## https://stat.ethz.ch/R-manual/R-devel/library/base/html/unique.html
probeIDs2GeneNames <- function(df){
  probeIDs <- df$Predictor
  # remove X affixes
  probeIDs <- sapply(na.omit(probeIDs),function(x) gsub("^X","", x))
  probeIDs <- sapply(na.omit(probeIDs),function(x) gsub("_max$","",x))
  geneNames <- toupper(unique(as.character(na.omit(annotate::lookUp(probeIDs, "hgu133a2.db", "SYMBOL")))))
  return(data.frame(Predictor = geneNames))
}

unigeneIDs2GeneNames <- function(df){
  unigeneIDs <- df$Predictor
  # if many reported, keep first
  unigeneIDs <- sapply(na.omit(unigeneIDs),function(x) gsub("^(Hs\\.[0-9]+).*","\\1",x))
  probeIDs <- unique(as.character(na.omit(AnnotationDbi::mget(unigeneIDs, ifnotfound=NA,envir=AnnotationDbi::revmap(hgu133a2UNIGENE)))))
  geneNames <- toupper(unique(as.character(na.omit(annotate::lookUp(probeIDs,"hgu133a2.db","SYMBOL")))))
  return(data.frame(Predictor = geneNames))
}


# gloabal variables
P_VAL_THRESH = 0.005

ENTREZ_UNIVERSE <- AnnotationDbi::select(hgu133a2.db,keys=ls(hgu133a2ENTREZID),columns="ENTREZID") %>% .[,"ENTREZID"]
ENTREZ_UNIVERSE <- ENTREZ_UNIVERSE[!duplicated(ENTREZ_UNIVERSE)]

DOWNLOAD_DIR <- "../data/downloads/"
RESULTS_DIR <- "../results/"
FIGURES_DIR <- paste(RESULTS_DIR,"Figs/",sep="")
GENELIST_DIR <- paste(RESULTS_DIR,"GeneLists/",sep="")

JH_S2_T0 <- paste(DOWNLOAD_DIR, "Subchallenge2_jhou_Time0_Predictors.csv",sep="")
JH_S2_T2 <- paste(DOWNLOAD_DIR, "Subchallenge2_jhou_Time24_Predictors.csv",sep="")
JH_S3_T0 <- paste(DOWNLOAD_DIR, "Subchallenge3_jhou_Time0_Predictors.csv",sep="")
JH_S3_T2 <- paste(DOWNLOAD_DIR, "Subchallenge3_jhou_Time24_Predictors.csv",sep="")

JH_S2_predictors <- read.csv(JH_S2_T0,header = TRUE,sep=",") %>%
  dplyr::transmute(Predictor = X) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()
JH_S3_predictors <- read.csv(JH_S3_T0,header = TRUE,sep=",") %>%
  dplyr::transmute(Predictor = X) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()

JH_S2_predictors <- read.csv(JH_S2_T2,header = TRUE,sep=",") %>%
  dplyr::transmute(Predictor = X) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(JH_S2_predictors)
JH_S3_predictors <- read.csv(JH_S3_T2,header = TRUE,sep=",") %>%
  dplyr::transmute(Predictor = X) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(JH_S3_predictors)

AY_S2_T0 <- paste(DOWNLOAD_DIR,"Subchallenge2_aydin_Time0_Predictors.csv",sep="")
AY_S2_T2 <- paste(DOWNLOAD_DIR,"Subchallenge2_aydin_Time24_Predictors.csv",sep="")
AY_S3_T0 <- paste(DOWNLOAD_DIR,"Subchallenge3_aydin_Time0_Predictors.csv",sep="")
AY_S3_T2 <- paste(DOWNLOAD_DIR,"Subchallenge3_aydin_Time24_Predictors.csv",sep="")

AY_S2_predictors <- read.csv(AY_S2_T0,sep=",") %>%
  dplyr::select(PROBE_SET_ID) %>%
  dplyr::transmute(Predictor = PROBE_SET_ID) %>%
  probeIDs2GeneNames()
AY_S3_predictors <- read.csv(AY_S3_T0,sep=",",header = FALSE) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames()

AY_S2_predictors <- read.csv(AY_S2_T2,sep=",") %>%
  dplyr::select(PROBE_SET_ID) %>%
  dplyr::transmute(Predictor = PROBE_SET_ID) %>%
  probeIDs2GeneNames() %>%
  rbind(AY_S2_predictors)
AY_S3_predictors <- read.csv(AY_S3_T2,sep=",",header = FALSE) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(AY_S3_predictors)

SS_S2_T0 <- paste(DOWNLOAD_DIR,"Subchallenge2_SSN DREAM TEAM-0_Predictors.csv",sep="")
SS_S2_T2 <- paste(DOWNLOAD_DIR,"Subchallenge2_SSN DREAM TEAM-24_Predictors.csv",sep="")
SS_S3_T0 <- paste(DOWNLOAD_DIR,"Subchallenge3_SSN DREAM TEAM-0_Predictors.csv",sep="")
SS_S3_T2 <- paste(DOWNLOAD_DIR,"Subchallenge3_SSN DREAM TEAM-24_Predictors.csv",sep="")

SS_S2_predictors <- read.csv(SS_S2_T0,sep=",",header=FALSE) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames()
SS_S3_predictors <- read.csv(SS_S3_T0,sep=",",header=FALSE) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames()

SS_S2_predictors <- read.csv(SS_S2_T2,sep=",",header=FALSE) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(SS_S2_predictors)
SS_S3_predictors <- read.csv(SS_S3_T2,sep=",",header=FALSE) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(SS_S3_predictors)

CW_S2_T2 <- paste(DOWNLOAD_DIR,"Subchallenge2_cwruPatho_Time24_Predictors.csv",sep="")
CW_S3_T0 <- paste(DOWNLOAD_DIR,"Subchallenge3_cwruPatho_Time0_Predictors.csv",sep="")
CW_S3_T2 <- paste(DOWNLOAD_DIR,"Subchallenge3_cwruPatho_Time24_Predictors.csv",sep="")

CW_S2_predictors <- read.csv(CW_S2_T2,sep=",") %>%
  dplyr::select(Probe.Set.ID) %>%
  dplyr::transmute(Predictor = Probe.Set.ID) %>%
  probeIDs2GeneNames()

CW_S3_predictors <- read.csv(CW_S3_T0,sep=",") %>%
  dplyr::select(Probe.Set.ID) %>%
  dplyr::transmute(Predictor = Probe.Set.ID) %>%
  probeIDs2GeneNames()
CW_S3_predictors <- read.csv(CW_S3_T2,sep=",") %>%
  dplyr::select(Probe.Set.ID) %>%
  dplyr::transmute(Predictor = Probe.Set.ID) %>%
  probeIDs2GeneNames() %>%
  rbind(CW_S3_predictors)

# calculate support

AY_S2_predictors <- AY_S2_predictors[!duplicated(AY_S2_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

SS_S2_predictors <- SS_S2_predictors[!duplicated(SS_S2_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

CW_S2_predictors <- CW_S2_predictors[!duplicated(CW_S2_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

JH_S2_predictors <- JH_S2_predictors[!duplicated(JH_S2_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

list.gene.sets <- c(
                    AY_S2_predictors,
                    SS_S2_predictors,
                    CW_S2_predictors,
                    JH_S2_predictors
                    )

names(list.gene.sets) <- c(
                           "AY",
                           "SS",
                           "CW",
                           "JH"
                           )

# from https://cran.r-project.org/web/packages/SuperExactTest/vignettes/set_html.html
length.gene.sets <- sapply(list.gene.sets,
                           FUN = length)

total=ENTREZ_UNIVERSE %>% length() # 24515

num.expected.overlap=total*do.call(base::prod,as.list(length.gene.sets/total))

sapply(0:min(length.gene.sets),function(i) SuperExactTest::dpsets(i, length.gene.sets, n=total))

res=SuperExactTest::supertest(list.gene.sets, n=total)

svg(filename=paste(FIGURES_DIR,"S2_Model_Intersections_(Spiral).svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
plot(res, sort.by="p-value")
dev.off()

svg(filename=paste(FIGURES_DIR,"S2_Model_Intersections_(Bar_Graph).svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
plot(res, Layout="landscape",degree=2:4, sort.by="p-value")
dev.off()

sigResDF <- summary(res)[["Table"]] %>%
  as.data.frame() %>%
  dplyr::filter(P.value < P_VAL_THRESH)

# end from https://cran.r-project.org/web/packages/SuperExactTest/vignettes/set_html.html

S2_predictors <- sigResDF %>%
  dplyr::select(Elements) %>%
  .[,] %>%
  strsplit(.,split=", ") %>%
  unlist() %>%
  unique()

AY_S3_predictors <- AY_S3_predictors[!duplicated(AY_S3_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

SS_S3_predictors <- SS_S3_predictors[!duplicated(SS_S3_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

CW_S3_predictors <- CW_S3_predictors[!duplicated(CW_S3_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

JH_S3_predictors <- JH_S3_predictors[!duplicated(JH_S3_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

list.gene.sets <- c(
                    AY_S3_predictors,
                    SS_S3_predictors,
                    CW_S3_predictors,
                    JH_S3_predictors
                    )

names(list.gene.sets) <- c(
                           "AY",
                           "SS",
                           "CW",
                           "JH"
                           )

# from https://cran.r-project.org/web/packages/SuperExactTest/vignettes/set_html.html
length.gene.sets <- sapply(list.gene.sets,
                           FUN = length)

total=ENTREZ_UNIVERSE %>% length() # 24515

num.expected.overlap=total*do.call(base::prod,as.list(length.gene.sets/total))

sapply(0:min(length.gene.sets),function(i) SuperExactTest::dpsets(i, length.gene.sets, n=total))

res=SuperExactTest::supertest(list.gene.sets, n=total)

svg(filename=paste(FIGURES_DIR,"S3_Model_Intersections_(Spiral).svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
plot(res, sort.by="p-value")
dev.off()

svg(filename=paste(FIGURES_DIR,"S3_Model_Intersections_(Bar_Graph).svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
plot(res, Layout="landscape",degree=2:4, sort.by="p-value")
dev.off()

sigResDF <- summary(res)[["Table"]] %>%
  as.data.frame() %>%
  dplyr::filter(P.value < P_VAL_THRESH)

# end from https://cran.r-project.org/web/packages/SuperExactTest/vignettes/set_html.html

S3_predictors <- sigResDF %>%
  dplyr::select(Elements) %>%
  .[,] %>%
  strsplit(.,split=", ") %>%
  unlist() %>%
  unique()

# reactome pathway enrichment analysis

## S2
title <- paste("S2 Predictors With p-value < ",P_VAL_THRESH,sep="")
S2_predictors %>%
  write(paste(GENELIST_DIR,"S2_Predictors.txt",sep=""))
svg(filename=paste(FIGURES_DIR,"S2_Predictor_PA.svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
S2_predictors %>%
  ReactomeBarplot(query_gene_symbols = .,
                  bar_color = "lightskyblue",
                  plot_title = title)
dev.off()

## S3
title <- paste("S3 Predictors With p-value < ",P_VAL_THRESH,sep="")
S3_predictors %>%
  write(paste(GENELIST_DIR,"S3_Predictors.txt",sep=""))
svg(filename=paste(FIGURES_DIR,"S3_Predictor_PA.svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
S3_predictors %>%
  ReactomeBarplot(query_gene_symbols = .,
                  bar_color = "midnightblue",
                  plot_title = title)
dev.off()

# Intersections

S2_S3_predictors <- SuperExactTest::intersect(S2_predictors,
                                                  S3_predictors)
S2_S3_predictors %>%
  write(paste(GENELIST_DIR,"S2_S3_Intersection.txt",sep=""))

S2_Only <- dplyr::setdiff(S2_predictors,S3_predictors)
S2_Only %>%
  write(paste(GENELIST_DIR,"S2_Only.txt",sep=""))

S3_Only <- dplyr::setdiff(S3_predictors,S2_predictors)
S3_Only %>%
  write(paste(GENELIST_DIR,"S3_Only.txt",sep=""))

svg(filename=paste(FIGURES_DIR,"PerSubchallengeVennDiagram.svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
VennDiagram::draw.pairwise.venn(area1=length(S2_predictors),
                              area2=length(S3_predictors),
                              category = c("S2","S3"),
                              cross.area=length(S2_S3_predictors),
                              fill=c("lightskyblue","midnightblue"))
dev.off()

## S2 S3
svg(filename=paste(FIGURES_DIR,"S2_S3_Predictor_PA.svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
S2_S3_predictors %>%
  ReactomeBarplot(query_gene_symbols = .,
                  bar_color = "blue",
                  plot_title = paste("S2 & S3 Predictors With p-value < ",P_VAL_THRESH))
dev.off()


## S2 only
svg(filename=paste(FIGURES_DIR,"S2-Only_Predictor_PA.svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
S2_Only %>%
  ReactomeBarplot(query_gene_symbols = .,
                  bar_color = "lightskyblue",
                  plot_title = paste("S2-Only Predictors with p-value < ",P_VAL_THRESH))
dev.off()

## S3 only
svg(filename=paste(FIGURES_DIR,"S3-Only_Predictor_PA.svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
S3_Only %>%
  ReactomeBarplot(query_gene_symbols = .,
                  bar_color = "midnightblue",
                  plot_title = paste("S3-Only Predictors with p-value < ",P_VAL_THRESH))
dev.off()

system_command = paste("qlmanage -t -s 1000 -o ", FIGURES_DIR, " ",FIGURES_DIR,"*.svg",sep="")
system(paste(system_command))
