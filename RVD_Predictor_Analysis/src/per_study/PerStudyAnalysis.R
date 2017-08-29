# ---
# title: Per Study Analysis
# author: Joshua Burkhart
# date: Aug 10, 2017
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

# gloabal variables
P_VAL_THRESH = 0.005

ENTREZ_UNIVERSE <- AnnotationDbi::select(hgu133a2.db,keys=ls(hgu133a2ENTREZID),columns="ENTREZID") %>% .[,"ENTREZID"]
ENTREZ_UNIVERSE <- ENTREZ_UNIVERSE[!duplicated(ENTREZ_UNIVERSE)]

DOWNLOAD_DIR <- "../data/downloads/"
RESULTS_DIR <- "../results/"
FIGURES_DIR <- paste(RESULTS_DIR,"Figs/",sep="")
GENELIST_DIR <- paste(RESULTS_DIR,"GeneLists/",sep="")

NA_S1_T0 <- paste(DOWNLOAD_DIR,"Subchallenge1_Nautilus_Time0_Predictors.csv",sep="")
NA_S1_T2 <- paste(DOWNLOAD_DIR,"Subchallenge1_Nautilus_Time24_Predictors.csv",sep="")
NA_S2_T0 <- paste(DOWNLOAD_DIR,"Subchallenge2_Nautilus_Time0_Predictors.csv",sep="")
NA_S2_T2 <- paste(DOWNLOAD_DIR,"Subchallenge2_Nautilus_Time24_Predictors.csv",sep="")

AY_S1_T0 <- paste(DOWNLOAD_DIR,"Subchallenge1_aydin_Time0_Predictors.csv",sep="")
AY_S1_T2 <- paste(DOWNLOAD_DIR,"Subchallenge1_aydin_Time24_Predictors.csv",sep="")
AY_S2_T0 <- paste(DOWNLOAD_DIR,"Subchallenge2_aydin_Time0_Predictors.csv",sep="")
AY_S2_T2 <- paste(DOWNLOAD_DIR,"Subchallenge2_aydin_Time24_Predictors.csv",sep="")
AY_S3_T0 <- paste(DOWNLOAD_DIR,"Subchallenge3_aydin_Time0_Predictors.csv",sep="")
AY_S3_T2 <- paste(DOWNLOAD_DIR,"Subchallenge3_aydin_Time24_Predictors.csv",sep="")

SS_S1_T0 <- paste(DOWNLOAD_DIR,"Subchallenge1_SSN DREAM TEAM_Time0_Predictors.csv",sep="")
SS_S1_T2 <- paste(DOWNLOAD_DIR,"Subchallenge1_SSN DREAM TEAM-24_Predictors.csv",sep="")
SS_S2_T0 <- paste(DOWNLOAD_DIR,"Subchallenge2_SSN DREAM TEAM-0_Predictors.csv",sep="")
SS_S2_T2 <- paste(DOWNLOAD_DIR,"Subchallenge2_SSN DREAM TEAM-24_Predictors.csv",sep="")
SS_S3_T0 <- paste(DOWNLOAD_DIR,"Subchallenge3_SSN DREAM TEAM-0_Predictors.csv",sep="")
SS_S3_T2 <- paste(DOWNLOAD_DIR,"Subchallenge3_SSN DREAM TEAM-24_Predictors.csv",sep="")

# Unigene IDs
# DEE4 = H1N1, DEE5 = H3N2, DUKE = Rhinovirus
TX_S1_T0 <- paste(DOWNLOAD_DIR,"SubChallenge1_TXsolo_Time0_Predictors.csv",sep="")
TX_S1_T2 <- paste(DOWNLOAD_DIR,"SubChallenge1_TXsolo_Time24_Predictors.csv",sep="")
TX_S2_T0 <- paste(DOWNLOAD_DIR,"SubChallenge2_TXsolo_Time0_Predictors.csv",sep="")
TX_S2_T2 <- paste(DOWNLOAD_DIR,"SubChallenge2_TXsolo_Time24_Predictors.csv",sep="")

CW_S1_T0 <- paste(DOWNLOAD_DIR,"Subchallenge1_cwruPatho_Time0_Predictors.csv",sep="")
CW_S1_T2 <- paste(DOWNLOAD_DIR,"Subchallenge1_cwruPatho_Time24_Predictors.csv",sep="")
CW_S2_T2 <- paste(DOWNLOAD_DIR,"Subchallenge2_cwruPatho_Time24_Predictors.csv",sep="")
CW_S3_T0 <- paste(DOWNLOAD_DIR,"Subchallenge3_cwruPatho_Time0_Predictors.csv",sep="")
CW_S3_T2 <- paste(DOWNLOAD_DIR,"Subchallenge3_cwruPatho_Time24_Predictors.csv",sep="")

#HGNC IDs
AG_S1_T0 <- paste(DOWNLOAD_DIR,"Subchallenge1_Aganita_Time0_Predictors.csv",sep="")
AG_S1_T2 <- paste(DOWNLOAD_DIR,"Subchallenge1_Aganita_Time24_Predictors.csv",sep="")
AG_S2_T0 <- paste(DOWNLOAD_DIR,"Subchallenge2_Aganita_Time0_Predictors.csv",sep="")
AG_S2_T2 <- paste(DOWNLOAD_DIR,"Subchallenge2_Aganita_Time24_Predictors.csv",sep="")
AG_S3_T0 <- paste(DOWNLOAD_DIR,"Subchallenge3_Aganita_Time0_Predictors.csv",sep="")
AG_S3_T2 <- paste(DOWNLOAD_DIR,"Subchallenge3_Aganita_Time24_Predictors.csv",sep="")

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


# H1N1

## Nautilus
NA_H1N1_predictors <- read.csv(NA_S1_T0,sep=";") %>%
  dplyr::filter(grepl("H1N1",Ensemble)) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()
NA_H1N1_predictors <- read.csv(NA_S1_T2,sep=";") %>%
  dplyr::filter(grepl("H1N1",Ensemble)) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(NA_H1N1_predictors)
NA_H1N1_predictors <- read.csv(NA_S2_T0,sep=";") %>%
  dplyr::filter(grepl("H1N1",Ensemble)) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(NA_H1N1_predictors)
NA_H1N1_predictors <- read.csv(NA_S2_T2,sep=";") %>%
  dplyr::filter(grepl("H1N1",Ensemble)) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(NA_H1N1_predictors)

## aydin
AY_H1N1_predictors <- read.csv(AY_S1_T0,sep=",") %>%
  dplyr::filter(grepl("H1N1",VIRUS_TYPE)) %>%
  dplyr::select(PROBE_SET_ID) %>%
  dplyr::transmute(Predictor = PROBE_SET_ID) %>%
  probeIDs2GeneNames()
AY_H1N1_predictors <- read.csv(AY_S1_T2,sep=",") %>%
  dplyr::filter(grepl("H1N1",VIRUS_TYPE)) %>%
  dplyr::select(PROBE_SET_ID) %>%
  dplyr::transmute(Predictor = PROBE_SET_ID) %>%
  probeIDs2GeneNames() %>%
  rbind(AY_H1N1_predictors)
AY_H1N1_predictors <- read.csv(AY_S2_T0,sep=",") %>%
  dplyr::filter(grepl("H1N1",VIRUS_TYPE)) %>%
  dplyr::select(PROBE_SET_ID) %>%
  dplyr::transmute(Predictor = PROBE_SET_ID) %>%
  probeIDs2GeneNames() %>%
  rbind(AY_H1N1_predictors)
AY_H1N1_predictors <- read.csv(AY_S2_T2,sep=",") %>%
  dplyr::filter(grepl("H1N1",VIRUS_TYPE)) %>%
  dplyr::select(PROBE_SET_ID) %>%
  dplyr::transmute(Predictor = PROBE_SET_ID) %>%
  probeIDs2GeneNames() %>%
  rbind(AY_H1N1_predictors)
AY_H1N1_predictors <- read.csv(AY_S3_T0,sep=",",header = FALSE) %>%
  dplyr::filter(grepl("H1N1",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(AY_H1N1_predictors)
AY_H1N1_predictors <- read.csv(AY_S3_T2,sep=",",header = FALSE) %>%
  dplyr::filter(grepl("H1N1",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(AY_H1N1_predictors)

## SSN_Dream_Team
SS_H1N1_predictors <- read.csv(SS_S1_T0,sep=",",header=FALSE) %>%
  dplyr::filter(grepl("H1N1",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames()
SS_H1N1_predictors <- read.csv(SS_S1_T2,sep=",",header=FALSE) %>%
  dplyr::filter(grepl("H1N1",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(SS_H1N1_predictors)
SS_H1N1_predictors <- read.csv(SS_S2_T0,sep=",",header=FALSE) %>%
  dplyr::filter(grepl("H1N1",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(SS_H1N1_predictors)
SS_H1N1_predictors <- read.csv(SS_S2_T2,sep=",",header=FALSE) %>%
  dplyr::filter(grepl("H1N1",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(SS_H1N1_predictors)
SS_H1N1_predictors <- read.csv(SS_S3_T0,sep=",",header=FALSE) %>%
  dplyr::filter(grepl("H1N1",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(SS_H1N1_predictors)
SS_H1N1_predictors <- read.csv(SS_S3_T2,sep=",",header=FALSE) %>%
  dplyr::filter(grepl("H1N1",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(SS_H1N1_predictors)

## Txsolo
TX_H1N1_predictors <- read.csv(TX_S1_T0,sep=",") %>%
  dplyr::filter(grepl("DEE4",virus)) %>%
  dplyr::select(feature) %>%
  dplyr::transmute(Predictor = feature) %>%
  unigeneIDs2GeneNames()
TX_H1N1_predictors <- read.csv(TX_S1_T2,sep=",") %>%
  dplyr::filter(grepl("DEE4",virus)) %>%
  dplyr::select(feature) %>%
  dplyr::transmute(Predictor = feature) %>%
  unigeneIDs2GeneNames() %>%
  rbind(TX_H1N1_predictors)
TX_H1N1_predictors <- read.csv(TX_S2_T0,sep=",") %>%
  dplyr::filter(grepl("DEE4",virus)) %>%
  dplyr::select(feature) %>%
  dplyr::transmute(Predictor = feature) %>%
  unigeneIDs2GeneNames() %>%
  rbind(TX_H1N1_predictors)
TX_H1N1_predictors <- read.csv(TX_S2_T2,sep=",") %>%
  dplyr::filter(grepl("DEE4",virus)) %>%
  dplyr::select(feature) %>%
  dplyr::transmute(Predictor = feature) %>%
  unigeneIDs2GeneNames() %>%
  rbind(TX_H1N1_predictors)

## cwruPatho
CW_H1N1_predictors <- read.csv(CW_S1_T0,sep=",") %>%
  dplyr::filter(grepl("H1N1",Virus)) %>%
  dplyr::select(Probe.Set.ID) %>%
  dplyr::transmute(Predictor = Probe.Set.ID) %>%
  probeIDs2GeneNames()
CW_H1N1_predictors <- read.csv(CW_S1_T2,sep=",") %>%
  dplyr::filter(grepl("H1N1",Virus)) %>%
  dplyr::select(Probe.Set.ID) %>%
  dplyr::transmute(Predictor = Probe.Set.ID) %>%
  probeIDs2GeneNames() %>%
  rbind(CW_H1N1_predictors)
CW_H1N1_predictors <- read.csv(CW_S2_T2,sep=",") %>%
  dplyr::filter(grepl("H1N1",Virus)) %>%
  dplyr::select(Probe.Set.ID) %>%
  dplyr::transmute(Predictor = Probe.Set.ID) %>%
  probeIDs2GeneNames() %>%
  rbind(CW_H1N1_predictors)
CW_H1N1_predictors <- read.csv(CW_S3_T0,sep=",") %>%
  dplyr::filter(grepl("H1N1",Virus)) %>%
  dplyr::select(Probe.Set.ID) %>%
  dplyr::transmute(Predictor = Probe.Set.ID) %>%
  probeIDs2GeneNames() %>%
  rbind(CW_H1N1_predictors)
CW_H1N1_predictors <- read.csv(CW_S3_T2,sep=",") %>%
  dplyr::filter(grepl("H1N1",Virus)) %>%
  dplyr::select(Probe.Set.ID) %>%
  dplyr::transmute(Predictor = Probe.Set.ID) %>%
  probeIDs2GeneNames() %>%
  rbind(CW_H1N1_predictors)

## Aganita
AG_H1N1_predictors <- read.csv(AG_S1_T0,sep=",") %>%
  dplyr::filter(grepl("H1N1",Study)) %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene)
AG_H1N1_predictors <- read.csv(AG_S1_T2,sep=",") %>%
  dplyr::filter(grepl("H1N1",Study)) %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene) %>%
  rbind(AG_H1N1_predictors)
AG_H1N1_predictors <- read.csv(AG_S2_T0,sep=",") %>%
  dplyr::filter(grepl("H1N1",Study)) %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene) %>%
  rbind(AG_H1N1_predictors)
AG_H1N1_predictors <- read.csv(AG_S2_T2,sep=",") %>%
  dplyr::filter(grepl("H1N1",Study)) %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene) %>%
  rbind(AG_H1N1_predictors)
AG_H1N1_predictors <- read.csv(AG_S3_T0,sep=",") %>%
  dplyr::filter(grepl("H1N1",Study)) %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene) %>%
  rbind(AG_H1N1_predictors)
AG_H1N1_predictors <- read.csv(AG_S3_T2,sep=",") %>%
  dplyr::filter(grepl("H1N1",Study)) %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene) %>%
  rbind(AG_H1N1_predictors)

# H3N2

## Nautilus
NA_H3N2_predictors <- read.csv(NA_S1_T0,sep=";") %>%
  dplyr::filter(grepl("H3N2",Ensemble)) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()
NA_H3N2_predictors <- read.csv(NA_S1_T2,sep=";") %>%
  dplyr::filter(grepl("H3N2",Ensemble)) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(NA_H3N2_predictors)
NA_H3N2_predictors <- read.csv(NA_S2_T0,sep=";") %>%
  dplyr::filter(grepl("H3N2",Ensemble)) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(NA_H3N2_predictors)
NA_H3N2_predictors <- read.csv(NA_S2_T2,sep=";") %>%
  dplyr::filter(grepl("H3N2",Ensemble)) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(NA_H3N2_predictors)

## aydin
AY_H3N2_predictors <- read.csv(AY_S1_T0,sep=",") %>%
  dplyr::filter(grepl("H3N2",VIRUS_TYPE)) %>%
  dplyr::select(PROBE_SET_ID) %>%
  dplyr::transmute(Predictor = PROBE_SET_ID) %>%
  probeIDs2GeneNames()
AY_H3N2_predictors <- read.csv(AY_S1_T2,sep=",") %>%
  dplyr::filter(grepl("H3N2",VIRUS_TYPE)) %>%
  dplyr::select(PROBE_SET_ID) %>%
  dplyr::transmute(Predictor = PROBE_SET_ID) %>%
  probeIDs2GeneNames() %>%
  rbind(AY_H3N2_predictors)
AY_H3N2_predictors <- read.csv(AY_S2_T0,sep=",") %>%
  dplyr::filter(grepl("H3N2",VIRUS_TYPE)) %>%
  dplyr::select(PROBE_SET_ID) %>%
  dplyr::transmute(Predictor = PROBE_SET_ID) %>%
  probeIDs2GeneNames() %>%
  rbind(AY_H3N2_predictors)
AY_H3N2_predictors <- read.csv(AY_S2_T2,sep=",") %>%
  dplyr::filter(grepl("H3N2",VIRUS_TYPE)) %>%
  dplyr::select(PROBE_SET_ID) %>%
  dplyr::transmute(Predictor = PROBE_SET_ID) %>%
  probeIDs2GeneNames() %>%
  rbind(AY_H3N2_predictors)
AY_H3N2_predictors <- read.csv(AY_S3_T0,sep=",",header = FALSE) %>%
  dplyr::filter(grepl("H3N2",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(AY_H3N2_predictors)
AY_H3N2_predictors <- read.csv(AY_S3_T2,sep=",",header = FALSE) %>%
  dplyr::filter(grepl("H3N2",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(AY_H3N2_predictors)

## SSN_Dream_Team
SS_H3N2_predictors <- read.csv(SS_S1_T0,sep=",",header=FALSE) %>%
  dplyr::filter(grepl("H3N2",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames()
SS_H3N2_predictors <- read.csv(SS_S1_T2,sep=",",header=FALSE) %>%
  dplyr::filter(grepl("H3N2",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(SS_H3N2_predictors)
SS_H3N2_predictors <- read.csv(SS_S2_T0,sep=",",header=FALSE) %>%
  dplyr::filter(grepl("H3N2",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(SS_H3N2_predictors)
SS_H3N2_predictors <- read.csv(SS_S2_T2,sep=",",header=FALSE) %>%
  dplyr::filter(grepl("H3N2",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(SS_H3N2_predictors)
SS_H3N2_predictors <- read.csv(SS_S3_T0,sep=",",header=FALSE) %>%
  dplyr::filter(grepl("H3N2",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(SS_H3N2_predictors)
SS_H3N2_predictors <- read.csv(SS_S3_T2,sep=",",header=FALSE) %>%
  dplyr::filter(grepl("H3N2",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(SS_H3N2_predictors)

## Txsolo
TX_H3N2_predictors <- read.csv(TX_S1_T0,sep=",") %>%
  dplyr::filter(grepl("DEE5",virus)) %>%
  dplyr::select(feature) %>%
  dplyr::transmute(Predictor = feature) %>%
  unigeneIDs2GeneNames()
TX_H3N2_predictors <- read.csv(TX_S1_T2,sep=",") %>%
  dplyr::filter(grepl("DEE5",virus)) %>%
  dplyr::select(feature) %>%
  dplyr::transmute(Predictor = feature) %>%
  unigeneIDs2GeneNames() %>%
  rbind(TX_H3N2_predictors)
TX_H3N2_predictors <- read.csv(TX_S2_T0,sep=",") %>%
  dplyr::filter(grepl("DEE5",virus)) %>%
  dplyr::select(feature) %>%
  dplyr::transmute(Predictor = feature) %>%
  unigeneIDs2GeneNames() %>%
  rbind(TX_H3N2_predictors)
TX_H3N2_predictors <- read.csv(TX_S2_T2,sep=",") %>%
  dplyr::filter(grepl("DEE5",virus)) %>%
  dplyr::select(feature) %>%
  dplyr::transmute(Predictor = feature) %>%
  unigeneIDs2GeneNames() %>%
  rbind(TX_H3N2_predictors)

## cwruPatho
CW_H3N2_predictors <- read.csv(CW_S1_T0,sep=",") %>%
  dplyr::filter(grepl("H3N2",Virus)) %>%
  dplyr::select(Probe.Set.ID) %>%
  dplyr::transmute(Predictor = Probe.Set.ID) %>%
  probeIDs2GeneNames()
CW_H3N2_predictors <- read.csv(CW_S1_T2,sep=",") %>%
  dplyr::filter(grepl("H3N2",Virus)) %>%
  dplyr::select(Probe.Set.ID) %>%
  dplyr::transmute(Predictor = Probe.Set.ID) %>%
  probeIDs2GeneNames() %>%
  rbind(CW_H3N2_predictors)
CW_H3N2_predictors <- read.csv(CW_S2_T2,sep=",") %>%
  dplyr::filter(grepl("H3N2",Virus)) %>%
  dplyr::select(Probe.Set.ID) %>%
  dplyr::transmute(Predictor = Probe.Set.ID) %>%
  probeIDs2GeneNames() %>%
  rbind(CW_H3N2_predictors)
CW_H3N2_predictors <- read.csv(CW_S3_T0,sep=",") %>%
  dplyr::filter(grepl("H3N2",Virus)) %>%
  dplyr::select(Probe.Set.ID) %>%
  dplyr::transmute(Predictor = Probe.Set.ID) %>%
  probeIDs2GeneNames() %>%
  rbind(CW_H3N2_predictors)
CW_H3N2_predictors <- read.csv(CW_S3_T2,sep=",") %>%
  dplyr::filter(grepl("H3N2",Virus)) %>%
  dplyr::select(Probe.Set.ID) %>%
  dplyr::transmute(Predictor = Probe.Set.ID) %>%
  probeIDs2GeneNames() %>%
  rbind(CW_H3N2_predictors)

## Aganita
AG_H3N2_predictors <- read.csv(AG_S1_T0,sep=",") %>%
  dplyr::filter(grepl("H3N2",Study)) %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene)
AG_H3N2_predictors <- read.csv(AG_S1_T2,sep=",") %>%
  dplyr::filter(grepl("H3N2",Study)) %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene) %>%
  rbind(AG_H3N2_predictors)
AG_H3N2_predictors <- read.csv(AG_S2_T0,sep=",") %>%
  dplyr::filter(grepl("H3N2",Study)) %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene) %>%
  rbind(AG_H3N2_predictors)
AG_H3N2_predictors <- read.csv(AG_S2_T2,sep=",") %>%
  dplyr::filter(grepl("H3N2",Study)) %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene) %>%
  rbind(AG_H3N2_predictors)
AG_H3N2_predictors <- read.csv(AG_S3_T0,sep=",") %>%
  dplyr::filter(grepl("H3N2",Study)) %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene) %>%
  rbind(AG_H3N2_predictors)
AG_H3N2_predictors <- read.csv(AG_S3_T2,sep=",") %>%
  dplyr::filter(grepl("H3N2",Study)) %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene) %>%
  rbind(AG_H3N2_predictors)

# Rhinovirus

## Nautilus
NA_Rhinovirus_predictors <- read.csv(NA_S1_T0,sep=";") %>%
  dplyr::filter(grepl("Rhinovirus",Ensemble)) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()
NA_Rhinovirus_predictors <- read.csv(NA_S1_T2,sep=";") %>%
  dplyr::filter(grepl("Rhinovirus",Ensemble)) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(NA_Rhinovirus_predictors)
NA_Rhinovirus_predictors <- read.csv(NA_S2_T0,sep=";") %>%
  dplyr::filter(grepl("Rhinovirus",Ensemble)) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(NA_Rhinovirus_predictors)
NA_Rhinovirus_predictors <- read.csv(NA_S2_T2,sep=";") %>%
  dplyr::filter(grepl("Rhinovirus",Ensemble)) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(NA_Rhinovirus_predictors)

## aydin
AY_Rhinovirus_predictors <- read.csv(AY_S1_T0,sep=",") %>%
  dplyr::filter(grepl("Rhinovirus",VIRUS_TYPE)) %>%
  dplyr::select(PROBE_SET_ID) %>%
  dplyr::transmute(Predictor = PROBE_SET_ID) %>%
  probeIDs2GeneNames()
AY_Rhinovirus_predictors <- read.csv(AY_S1_T2,sep=",") %>%
  dplyr::filter(grepl("Rhinovirus",VIRUS_TYPE)) %>%
  dplyr::select(PROBE_SET_ID) %>%
  dplyr::transmute(Predictor = PROBE_SET_ID) %>%
  probeIDs2GeneNames() %>%
  rbind(AY_Rhinovirus_predictors)
AY_Rhinovirus_predictors <- read.csv(AY_S2_T0,sep=",") %>%
  dplyr::filter(grepl("Rhinovirus",VIRUS_TYPE)) %>%
  dplyr::select(PROBE_SET_ID) %>%
  dplyr::transmute(Predictor = PROBE_SET_ID) %>%
  probeIDs2GeneNames() %>%
  rbind(AY_Rhinovirus_predictors)
AY_Rhinovirus_predictors <- read.csv(AY_S2_T2,sep=",") %>%
  dplyr::filter(grepl("Rhinovirus",VIRUS_TYPE)) %>%
  dplyr::select(PROBE_SET_ID) %>%
  dplyr::transmute(Predictor = PROBE_SET_ID) %>%
  probeIDs2GeneNames() %>%
  rbind(AY_Rhinovirus_predictors)
AY_Rhinovirus_predictors <- read.csv(AY_S3_T0,sep=",",header = FALSE) %>%
  dplyr::filter(grepl("Rhinovirus",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(AY_Rhinovirus_predictors)
AY_Rhinovirus_predictors <- read.csv(AY_S3_T2,sep=",",header = FALSE) %>%
  dplyr::filter(grepl("Rhinovirus",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(AY_Rhinovirus_predictors)

## SSN_Dream_Team
SS_Rhinovirus_predictors <- read.csv(SS_S1_T0,sep=",",header=FALSE) %>%
  dplyr::filter(grepl("Rhinovirus",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames()
SS_Rhinovirus_predictors <- read.csv(SS_S1_T2,sep=",",header=FALSE) %>%
  dplyr::filter(grepl("Rhinovirus",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(SS_Rhinovirus_predictors)
SS_Rhinovirus_predictors <- read.csv(SS_S2_T0,sep=",",header=FALSE) %>%
  dplyr::filter(grepl("Rhinovirus",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(SS_Rhinovirus_predictors)
SS_Rhinovirus_predictors <- read.csv(SS_S2_T2,sep=",",header=FALSE) %>%
  dplyr::filter(grepl("Rhinovirus",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(SS_Rhinovirus_predictors)
SS_Rhinovirus_predictors <- read.csv(SS_S3_T0,sep=",",header=FALSE) %>%
  dplyr::filter(grepl("Rhinovirus",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(SS_Rhinovirus_predictors)
SS_Rhinovirus_predictors <- read.csv(SS_S3_T2,sep=",",header=FALSE) %>%
  dplyr::filter(grepl("Rhinovirus",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(SS_Rhinovirus_predictors)

## Txsolo
TX_Rhinovirus_predictors <- read.csv(TX_S1_T0,sep=",") %>%
  dplyr::filter(grepl("DUKE",virus)) %>%
  dplyr::select(feature) %>%
  dplyr::transmute(Predictor = feature) %>%
  unigeneIDs2GeneNames()
TX_Rhinovirus_predictors <- read.csv(TX_S1_T2,sep=",") %>%
  dplyr::filter(grepl("DUKE",virus)) %>%
  dplyr::select(feature) %>%
  dplyr::transmute(Predictor = feature) %>%
  unigeneIDs2GeneNames() %>%
  rbind(TX_Rhinovirus_predictors)
TX_Rhinovirus_predictors <- read.csv(TX_S2_T0,sep=",") %>%
  dplyr::filter(grepl("DUKE",virus)) %>%
  dplyr::select(feature) %>%
  dplyr::transmute(Predictor = feature) %>%
  unigeneIDs2GeneNames() %>%
  rbind(TX_Rhinovirus_predictors)
TX_Rhinovirus_predictors <- read.csv(TX_S2_T2,sep=",") %>%
  dplyr::filter(grepl("DUKE",virus)) %>%
  dplyr::select(feature) %>%
  dplyr::transmute(Predictor = feature) %>%
  unigeneIDs2GeneNames() %>%
  rbind(TX_Rhinovirus_predictors)

## cwruPatho
CW_Rhinovirus_predictors <- read.csv(CW_S1_T0,sep=",") %>%
  dplyr::filter(grepl("Rhinovirus",Virus)) %>%
  dplyr::select(Probe.Set.ID) %>%
  dplyr::transmute(Predictor = Probe.Set.ID) %>%
  probeIDs2GeneNames()
CW_Rhinovirus_predictors <- read.csv(CW_S1_T2,sep=",") %>%
  dplyr::filter(grepl("Rhinovirus",Virus)) %>%
  dplyr::select(Probe.Set.ID) %>%
  dplyr::transmute(Predictor = Probe.Set.ID) %>%
  probeIDs2GeneNames() %>%
  rbind(CW_Rhinovirus_predictors)
CW_Rhinovirus_predictors <- read.csv(CW_S2_T2,sep=",") %>%
  dplyr::filter(grepl("Rhinovirus",Virus)) %>%
  dplyr::select(Probe.Set.ID) %>%
  dplyr::transmute(Predictor = Probe.Set.ID) %>%
  probeIDs2GeneNames() %>%
  rbind(CW_Rhinovirus_predictors)
CW_Rhinovirus_predictors <- read.csv(CW_S3_T0,sep=",") %>%
  dplyr::filter(grepl("Rhinovirus",Virus)) %>%
  dplyr::select(Probe.Set.ID) %>%
  dplyr::transmute(Predictor = Probe.Set.ID) %>%
  probeIDs2GeneNames() %>%
  rbind(CW_Rhinovirus_predictors)
CW_Rhinovirus_predictors <- read.csv(CW_S3_T2,sep=",") %>%
  dplyr::filter(grepl("Rhinovirus",Virus)) %>%
  dplyr::select(Probe.Set.ID) %>%
  dplyr::transmute(Predictor = Probe.Set.ID) %>%
  probeIDs2GeneNames() %>%
  rbind(CW_Rhinovirus_predictors)

## Aganita
AG_Rhinovirus_predictors <- read.csv(AG_S1_T0,sep=",") %>%
  dplyr::filter(grepl("Rhinovirus",Study)) %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene)
AG_Rhinovirus_predictors <- read.csv(AG_S1_T2,sep=",") %>%
  dplyr::filter(grepl("Rhinovirus",Study)) %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene) %>%
  rbind(AG_Rhinovirus_predictors)
AG_Rhinovirus_predictors <- read.csv(AG_S2_T0,sep=",") %>%
  dplyr::filter(grepl("Rhinovirus",Study)) %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene) %>%
  rbind(AG_Rhinovirus_predictors)
AG_Rhinovirus_predictors <- read.csv(AG_S2_T2,sep=",") %>%
  dplyr::filter(grepl("Rhinovirus",Study)) %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene) %>%
  rbind(AG_Rhinovirus_predictors)
AG_Rhinovirus_predictors <- read.csv(AG_S3_T0,sep=",") %>%
  dplyr::filter(grepl("Rhinovirus",Study)) %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene) %>%
  rbind(AG_Rhinovirus_predictors)
AG_Rhinovirus_predictors <- read.csv(AG_S3_T2,sep=",") %>%
  dplyr::filter(grepl("Rhinovirus",Study)) %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene) %>%
  rbind(AG_Rhinovirus_predictors)

# calculate support

## H1N1

NA_H1N1_predictors <- NA_H1N1_predictors[!duplicated(NA_H1N1_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")
AY_H1N1_predictors <- AY_H1N1_predictors[!duplicated(AY_H1N1_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")
SS_H1N1_predictors <- SS_H1N1_predictors[!duplicated(SS_H1N1_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")
TX_H1N1_predictors <- TX_H1N1_predictors[!duplicated(TX_H1N1_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")
CW_H1N1_predictors <- CW_H1N1_predictors[!duplicated(CW_H1N1_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")
AG_H1N1_predictors <- AG_H1N1_predictors[!duplicated(AG_H1N1_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

list.gene.sets <- c(NA_H1N1_predictors,
                    AY_H1N1_predictors,
                    SS_H1N1_predictors,
                    TX_H1N1_predictors,
                    CW_H1N1_predictors,
                    AG_H1N1_predictors)

names(list.gene.sets) <- c("NA","AY","SS","TX","CW","AG")

# from https://cran.r-project.org/web/packages/SuperExactTest/vignettes/set_html.html
length.gene.sets <- sapply(list.gene.sets,
                          FUN = length)

total=ENTREZ_UNIVERSE %>% length() # 24515

num.expected.overlap=total*do.call(base::prod,as.list(length.gene.sets/total))

sapply(0:min(length.gene.sets),function(i) SuperExactTest::dpsets(i, length.gene.sets, n=total))

res=SuperExactTest::supertest(list.gene.sets, n=total)

svg(filename=paste(FIGURES_DIR,"H1N1_Model_Intersections_(Spiral).svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
plot(res, sort.by="p-value")
dev.off()

svg(filename=paste(FIGURES_DIR,"H1N1_Model_Intersections_(Bar_Graph).svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
plot(res, Layout="landscape",degree=2:4, sort.by="p-value")
dev.off()

sigResDF <- summary(res)[["Table"]] %>%
  as.data.frame() %>%
  dplyr::filter(P.value < P_VAL_THRESH)

# end from https://cran.r-project.org/web/packages/SuperExactTest/vignettes/set_html.html

H1N1_predictors <- sigResDF %>%
  dplyr::select(Elements) %>%
  .[,] %>%
  strsplit(.,split=", ") %>%
  unlist() %>%
  unique()

## H3N2

NA_H3N2_predictors <- NA_H3N2_predictors[!duplicated(NA_H3N2_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")
AY_H3N2_predictors <- AY_H3N2_predictors[!duplicated(AY_H3N2_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")
SS_H3N2_predictors <- SS_H3N2_predictors[!duplicated(SS_H3N2_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")
TX_H3N2_predictors <- TX_H3N2_predictors[!duplicated(TX_H3N2_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")
CW_H3N2_predictors <- CW_H3N2_predictors[!duplicated(CW_H3N2_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")
AG_H3N2_predictors <- AG_H3N2_predictors[!duplicated(AG_H3N2_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

list.gene.sets <- c(NA_H3N2_predictors,
                    AY_H3N2_predictors,
                    SS_H3N2_predictors,
                    TX_H3N2_predictors,
                    CW_H3N2_predictors,
                    AG_H3N2_predictors)

names(list.gene.sets) <- c("NA","AY","SS","TX","CW","AG")

# from https://cran.r-project.org/web/packages/SuperExactTest/vignettes/set_html.html
length.gene.sets <- sapply(list.gene.sets,
                           FUN = length)

total=ENTREZ_UNIVERSE %>% length() # 24515

num.expected.overlap=total*do.call(base::prod,as.list(length.gene.sets/total))

sapply(0:min(length.gene.sets),function(i) SuperExactTest::dpsets(i, length.gene.sets, n=total))

res=SuperExactTest::supertest(list.gene.sets, n=total)

svg(filename=paste(FIGURES_DIR,"H3N2_Model_Intersections_(Spiral).svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
plot(res, sort.by="p-value")
dev.off()

svg(filename=paste(FIGURES_DIR,"H3N2_Model_Intersections_(Bar_Graph).svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
plot(res, Layout="landscape",degree=2:4, sort.by="p-value")
dev.off()

sigResDF <- summary(res)[["Table"]] %>%
  as.data.frame() %>%
  dplyr::filter(P.value < P_VAL_THRESH)

# end from https://cran.r-project.org/web/packages/SuperExactTest/vignettes/set_html.html

H3N2_predictors <- sigResDF %>%
  dplyr::select(Elements) %>%
  .[,] %>%
  strsplit(.,split=", ") %>%
  unlist() %>%
  unique()

## Rhinovirus

NA_Rhinovirus_predictors <- NA_Rhinovirus_predictors[!duplicated(NA_Rhinovirus_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")
AY_Rhinovirus_predictors <- AY_Rhinovirus_predictors[!duplicated(AY_Rhinovirus_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")
SS_Rhinovirus_predictors <- SS_Rhinovirus_predictors[!duplicated(SS_Rhinovirus_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")
TX_Rhinovirus_predictors <- TX_Rhinovirus_predictors[!duplicated(TX_Rhinovirus_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")
CW_Rhinovirus_predictors <- CW_Rhinovirus_predictors[!duplicated(CW_Rhinovirus_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")
AG_Rhinovirus_predictors <- AG_Rhinovirus_predictors[!duplicated(AG_Rhinovirus_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

list.gene.sets <- c(NA_Rhinovirus_predictors,
                    AY_Rhinovirus_predictors,
                    SS_Rhinovirus_predictors,
                    TX_Rhinovirus_predictors,
                    CW_Rhinovirus_predictors,
                    AG_Rhinovirus_predictors)

names(list.gene.sets) <- c("NA","AY","SS","TX","CW","AG")

# from https://cran.r-project.org/web/packages/SuperExactTest/vignettes/set_html.html
length.gene.sets <- sapply(list.gene.sets,
                           FUN = length)

total=ENTREZ_UNIVERSE %>% length() # 24515

num.expected.overlap=total*do.call(base::prod,as.list(length.gene.sets/total))

sapply(0:min(length.gene.sets),function(i) SuperExactTest::dpsets(i, length.gene.sets, n=total))

res=SuperExactTest::supertest(list.gene.sets, n=total)

svg(filename=paste(FIGURES_DIR,"Rhinovirus_Model_Intersections_(Spiral).svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
plot(res, sort.by="p-value")
dev.off()

svg(filename=paste(FIGURES_DIR,"Rhinovirus_Model_Intersections_(Bar_Graph).svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
plot(res, Layout="landscape",degree=2:4, sort.by="p-value")
dev.off()

sigResDF <- summary(res)[["Table"]] %>%
  as.data.frame() %>%
  dplyr::filter(P.value < P_VAL_THRESH)

# end from https://cran.r-project.org/web/packages/SuperExactTest/vignettes/set_html.html

Rhinovirus_predictors <- sigResDF %>%
  dplyr::select(Elements) %>%
  .[,] %>%
  strsplit(.,split=", ") %>%
  unlist() %>%
  unique()

# reactome pathway enrichment analysis

## H1N1
title <- paste("H1N1 Predictors With p-value < ",P_VAL_THRESH,sep="")
H1N1_predictors %>%
  write(paste(GENELIST_DIR,"H1N1_Predictors.txt",sep=""))
svg(filename=paste(FIGURES_DIR,"H1N1_Predictor_PA.svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
H1N1_predictors %>%
  ReactomeBarplot(query_gene_symbols = .,
                  bar_color = "pink",
                  plot_title = title)
dev.off()

## H3N2
title <- paste("H3N2 Predictors With p-value < ",P_VAL_THRESH,sep="")
H3N2_predictors %>%
  write(paste(GENELIST_DIR,"H3N2_Predictors.txt",sep=""))
svg(filename=paste(FIGURES_DIR,"H3N2_Predictor_PA.svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
H3N2_predictors %>%
  ReactomeBarplot(query_gene_symbols = .,
                  bar_color = "red",
                  plot_title = title)
dev.off()

## Rhinovirus
title <- paste("Rhinovirus Predictors With p-value < ",P_VAL_THRESH)
Rhinovirus_predictors %>%
  write(paste(GENELIST_DIR,"Rhinovirus_Predictors.txt",sep=""))
svg(filename=paste(FIGURES_DIR,"Rhinovirus_Predictor_PA.svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
Rhinovirus_predictors %>%
  ReactomeBarplot(query_gene_symbols = .,
                  bar_color = "gold",
                  plot_title = title)
dev.off()

# Intersections

H1N1_H3N2_predictors <- SuperExactTest::intersect(H1N1_predictors,
                                                  H3N2_predictors)
H1N1_H3N2_predictors %>%
  write(paste(GENELIST_DIR,"H1N1_H3N2_Intersection.txt",sep=""))
H3N2_Rhinovirus_predictors <- SuperExactTest::intersect(H3N2_predictors,
                                                        Rhinovirus_predictors)
H3N2_Rhinovirus_predictors %>%
  write(paste(GENELIST_DIR,"H3N2_Rhinovirus_Intersection.txt",sep=""))
H1N1_Rhinovirus_predictors <- SuperExactTest::intersect(H1N1_predictors,
                                                        Rhinovirus_predictors)
H1N1_Rhinovirus_predictors %>%
  write(paste(GENELIST_DIR,"H1N1_Rhinovirus_Intersection.txt",sep=""))
H1N1_H3N2_Rhinovirus_predictors <- SuperExactTest::intersect(H1N1_predictors,
                                                             H3N2_predictors,
                                                             Rhinovirus_predictors)
H1N1_H3N2_Rhinovirus_predictors %>%
  write(paste(GENELIST_DIR,"H1N1_H3N2_Rhinovirus_Intersection.txt",sep=""))

H1N1_Only <- dplyr::setdiff(H1N1_predictors,dplyr::union(H3N2_predictors,Rhinovirus_predictors))
H1N1_Only %>%
  write(paste(GENELIST_DIR,"H1N1_Only.txt",sep=""))

H3N2_Only <- dplyr::setdiff(H3N2_predictors,dplyr::union(H1N1_predictors,Rhinovirus_predictors))
H3N2_Only %>%
  write(paste(GENELIST_DIR,"H3N2_Only.txt",sep=""))

Rhinovirus_Only <- dplyr::setdiff(Rhinovirus_predictors,dplyr::union(H1N1_predictors,H3N2_predictors))
Rhinovirus_Only %>%
  write(paste(GENELIST_DIR,"Rhinovirus_Only.txt",sep=""))

svg(filename=paste(FIGURES_DIR,"PerStudyVennDiagram.svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
VennDiagram::draw.triple.venn(area1=length(H1N1_predictors),
                              area2=length(H3N2_predictors),
                              area3=length(Rhinovirus_predictors),
                              category = c("H1N1","H3N2","Rhinovirus"),
                              n12=length(H1N1_H3N2_predictors),
                              n23=length(H3N2_Rhinovirus_predictors),
                              n13=length(H1N1_Rhinovirus_predictors),
                              n123=length(H1N1_H3N2_Rhinovirus_predictors),
                              fill=c("pink","red","gold"))
dev.off()

## H1N1 H3N2
svg(filename=paste(FIGURES_DIR,"H1N1_H3N2_Predictor_PA.svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
H1N1_H3N2_predictors %>%
  ReactomeBarplot(query_gene_symbols = .,
                  bar_color = "indianred1",
                  plot_title = paste("H1N1 & H3N2 Predictors With p-value < ",P_VAL_THRESH))
dev.off()

## H3N2 Rhinovirus
svg(filename=paste(FIGURES_DIR,"H3N2_Rhinovirus_Predictor_PA.svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
H3N2_Rhinovirus_predictors %>%
  ReactomeBarplot(query_gene_symbols = .,
                  bar_color = "orange",
                  plot_title = paste("H3N2 & Rhinovirus Predictors With p-value < ",P_VAL_THRESH))
dev.off()

## H1N1 Rhinovirus
svg(filename=paste(FIGURES_DIR,"H1N1_Rhinovirus_Predictor_PA.svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
H1N1_Rhinovirus_predictors %>%
  ReactomeBarplot(query_gene_symbols = .,
                  bar_color = "sandybrown",
                  plot_title = paste("H1N1 & Rhinovirus Predictors With p-value < ",P_VAL_THRESH))
dev.off()

## H1N1 H3N2 Rhinovirus (should be GPX1)
svg(filename=paste(FIGURES_DIR,"H1N1_H3N2_Rhinovirus_Predictor_PA.svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
H1N1_H3N2_Rhinovirus_predictors %>% 
  ReactomeBarplot(query_gene_symbols = .,
                  bar_color = "lightsalmon",
                  plot_title = paste("H1N1 & H3N2 & Rhinovirus Predictors With p-value < ",P_VAL_THRESH))
dev.off()

## H1N1 only
svg(filename=paste(FIGURES_DIR,"H1N1-Only_Predictor_PA.svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
H1N1_Only %>%
  ReactomeBarplot(query_gene_symbols = .,
                  bar_color = "pink",
                  plot_title = paste("H1N1-Only Predictors with p-value < ",P_VAL_THRESH))
dev.off()

## H3N2 only
svg(filename=paste(FIGURES_DIR,"H3N2-Only_Predictor_PA.svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
H3N2_Only %>%
  ReactomeBarplot(query_gene_symbols = .,
                  bar_color = "red",
                  plot_title = paste("H3N2-Only Predictors with p-value < ",P_VAL_THRESH))
dev.off()

## Rhinovirus only
svg(filename=paste(FIGURES_DIR,"Rhinovirus-Only_Predictor_PA.svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
Rhinovirus_Only %>%
  ReactomeBarplot(query_gene_symbols = .,
                  bar_color = "gold",
                  plot_title = paste("Rhinovirus-Only Predictors with p-value < ",P_VAL_THRESH))
dev.off()

system_command = paste("qlmanage -t -s 1000 -o ", FIGURES_DIR, " ",FIGURES_DIR,"*.svg",sep="")
system(paste(system_command))
