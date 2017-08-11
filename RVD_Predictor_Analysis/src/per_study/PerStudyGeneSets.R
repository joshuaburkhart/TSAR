# ---
# title: Per Study Gene Sets
# author: Joshua Burkhart
# date: Aug 10, 2017
# ---

# load libraries
library(SuperExactTest)
library(org.Hs.eg.db)
library(hgu133a2.db)
library(ReactomePA)
library(annotate)
library(magrittr)
library(ggplot2)
library(dplyr)

# gloabal variables
ENTREZ_UNIVERSE <- AnnotationDbi::select(hgu133a2.db,keys=ls(hgu133a2ENTREZID),columns="ENTREZID") %>% .[,"ENTREZID"]

DOWNLOAD_DIR <- "/home/burkhart/Software/TSAR/RVD_Predictor_Analysis/data/downloads/"
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
  
  entrez_universe <- ENTREZ_UNIVERSE
  print(query_gene_symbols)
  entrez_query <- AnnotationDbi::select(org.Hs.eg.db, query_gene_symbols, "ENTREZID", "SYMBOL")
  
  entrez_universe <- entrez_universe[!duplicated(entrez_universe)]
  entrez_query <- entrez_query[!duplicated(entrez_query[,1]), 2]
  
  reactomeEnrich <- ReactomePA::enrichPathway(gene=unique(entrez_query),
                                  pvalueCutoff = 0.05,
                                  readable = T)
  reactomeEnrichAtLeast2GenesPerPathway <- reactomeEnrich %>% 
    as.data.frame() %>%
    dplyr::filter(Count >= 2)
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
    labs(x="Reactome Pathway",y="-log10(p-value)",title=plot_title) %>%
    return()
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
SS_H1N1_predictors <- read.csv(SS_S3_T0,sep=",",header=FALSE) %>%
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
H3N2_predictors <- read.csv(NA_S1_T0,sep=";") %>%
  dplyr::filter(grepl("H3N2",Ensemble)) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()
H3N2_predictors <- read.csv(NA_S1_T2,sep=";") %>%
  dplyr::filter(grepl("H3N2",Ensemble)) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(H3N2_predictors)
H3N2_predictors <- read.csv(NA_S2_T0,sep=";") %>%
  dplyr::filter(grepl("H3N2",Ensemble)) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(H3N2_predictors)
H3N2_predictors <- read.csv(NA_S2_T2,sep=";") %>%
  dplyr::filter(grepl("H3N2",Ensemble)) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(H3N2_predictors)

## aydin
H3N2_predictors <- read.csv(AY_S1_T0,sep=",") %>%
  dplyr::filter(grepl("H3N2",VIRUS_TYPE)) %>%
  dplyr::select(PROBE_SET_ID) %>%
  dplyr::transmute(Predictor = PROBE_SET_ID) %>%
  probeIDs2GeneNames() %>%
  rbind(H3N2_predictors)
H3N2_predictors <- read.csv(AY_S1_T2,sep=",") %>%
  dplyr::filter(grepl("H3N2",VIRUS_TYPE)) %>%
  dplyr::select(PROBE_SET_ID) %>%
  dplyr::transmute(Predictor = PROBE_SET_ID) %>%
  probeIDs2GeneNames() %>%
  rbind(H3N2_predictors)
H3N2_predictors <- read.csv(AY_S2_T0,sep=",") %>%
  dplyr::filter(grepl("H3N2",VIRUS_TYPE)) %>%
  dplyr::select(PROBE_SET_ID) %>%
  dplyr::transmute(Predictor = PROBE_SET_ID) %>%
  probeIDs2GeneNames() %>%
  rbind(H3N2_predictors)
H3N2_predictors <- read.csv(AY_S2_T2,sep=",") %>%
  dplyr::filter(grepl("H3N2",VIRUS_TYPE)) %>%
  dplyr::select(PROBE_SET_ID) %>%
  dplyr::transmute(Predictor = PROBE_SET_ID) %>%
  probeIDs2GeneNames() %>%
  rbind(H3N2_predictors)
H3N2_predictors <- read.csv(AY_S3_T0,sep=",",header = FALSE) %>%
  dplyr::filter(grepl("H3N2",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(H3N2_predictors)
H3N2_predictors <- read.csv(AY_S3_T2,sep=",",header = FALSE) %>%
  dplyr::filter(grepl("H3N2",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(H3N2_predictors)

## SSN_Dream_Team
H3N2_predictors <- read.csv(SS_S1_T0,sep=",",header=FALSE) %>%
  dplyr::filter(grepl("H3N2",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(H3N2_predictors)
H3N2_predictors <- read.csv(SS_S1_T2,sep=",",header=FALSE) %>%
  dplyr::filter(grepl("H3N2",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(H3N2_predictors)
H3N2_predictors <- read.csv(SS_S2_T0,sep=",",header=FALSE) %>%
  dplyr::filter(grepl("H3N2",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(H3N2_predictors)
H3N2_predictors <- read.csv(SS_S2_T2,sep=",",header=FALSE) %>%
  dplyr::filter(grepl("H3N2",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(H3N2_predictors)
H3N2_predictors <- read.csv(SS_S3_T0,sep=",",header=FALSE) %>%
  dplyr::filter(grepl("H3N2",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(H3N2_predictors)
H3N2_predictors <- read.csv(SS_S3_T0,sep=",",header=FALSE) %>%
  dplyr::filter(grepl("H3N2",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(H3N2_predictors)

## Txsolo
H3N2_predictors <- read.csv(TX_S1_T0,sep=",") %>%
  dplyr::filter(grepl("DEE5",virus)) %>%
  dplyr::select(feature) %>%
  dplyr::transmute(Predictor = feature) %>%
  unigeneIDs2GeneNames() %>%
  rbind(H3N2_predictors)
H3N2_predictors <- read.csv(TX_S1_T2,sep=",") %>%
  dplyr::filter(grepl("DEE5",virus)) %>%
  dplyr::select(feature) %>%
  dplyr::transmute(Predictor = feature) %>%
  unigeneIDs2GeneNames() %>%
  rbind(H3N2_predictors)
H3N2_predictors <- read.csv(TX_S2_T0,sep=",") %>%
  dplyr::filter(grepl("DEE5",virus)) %>%
  dplyr::select(feature) %>%
  dplyr::transmute(Predictor = feature) %>%
  unigeneIDs2GeneNames() %>%
  rbind(H3N2_predictors)
H3N2_predictors <- read.csv(TX_S2_T2,sep=",") %>%
  dplyr::filter(grepl("DEE5",virus)) %>%
  dplyr::select(feature) %>%
  dplyr::transmute(Predictor = feature) %>%
  unigeneIDs2GeneNames() %>%
  rbind(H3N2_predictors)

## cwruPatho
H3N2_predictors <- read.csv(CW_S1_T0,sep=",") %>%
  dplyr::filter(grepl("H3N2",Virus)) %>%
  dplyr::select(Probe.Set.ID) %>%
  dplyr::transmute(Predictor = Probe.Set.ID) %>%
  probeIDs2GeneNames() %>%
  rbind(H3N2_predictors)
H3N2_predictors <- read.csv(CW_S1_T2,sep=",") %>%
  dplyr::filter(grepl("H3N2",Virus)) %>%
  dplyr::select(Probe.Set.ID) %>%
  dplyr::transmute(Predictor = Probe.Set.ID) %>%
  probeIDs2GeneNames() %>%
  rbind(H3N2_predictors)
H3N2_predictors <- read.csv(CW_S2_T2,sep=",") %>%
  dplyr::filter(grepl("H3N2",Virus)) %>%
  dplyr::select(Probe.Set.ID) %>%
  dplyr::transmute(Predictor = Probe.Set.ID) %>%
  probeIDs2GeneNames() %>%
  rbind(H3N2_predictors)
H3N2_predictors <- read.csv(CW_S3_T0,sep=",") %>%
  dplyr::filter(grepl("H3N2",Virus)) %>%
  dplyr::select(Probe.Set.ID) %>%
  dplyr::transmute(Predictor = Probe.Set.ID) %>%
  probeIDs2GeneNames() %>%
  rbind(H3N2_predictors)
H3N2_predictors <- read.csv(CW_S3_T2,sep=",") %>%
  dplyr::filter(grepl("H3N2",Virus)) %>%
  dplyr::select(Probe.Set.ID) %>%
  dplyr::transmute(Predictor = Probe.Set.ID) %>%
  probeIDs2GeneNames() %>%
  rbind(H3N2_predictors)

## Aganita
H3N2_predictors <- read.csv(AG_S1_T0,sep=",") %>%
  dplyr::filter(grepl("H3N2",Study)) %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene) %>%
  rbind(H3N2_predictors)
H3N2_predictors <- read.csv(AG_S1_T2,sep=",") %>%
  dplyr::filter(grepl("H3N2",Study)) %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene) %>%
  rbind(H3N2_predictors)
H3N2_predictors <- read.csv(AG_S2_T0,sep=",") %>%
  dplyr::filter(grepl("H3N2",Study)) %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene) %>%
  rbind(H3N2_predictors)
H3N2_predictors <- read.csv(AG_S2_T2,sep=",") %>%
  dplyr::filter(grepl("H3N2",Study)) %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene) %>%
  rbind(H3N2_predictors)
H3N2_predictors <- read.csv(AG_S3_T0,sep=",") %>%
  dplyr::filter(grepl("H3N2",Study)) %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene) %>%
  rbind(H3N2_predictors)
H3N2_predictors <- read.csv(AG_S3_T2,sep=",") %>%
  dplyr::filter(grepl("H3N2",Study)) %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene) %>%
  rbind(H3N2_predictors)

# Rhinovirus

## Nautilus
Rhinovirus_predictors <- read.csv(NA_S1_T0,sep=";") %>%
  dplyr::filter(grepl("Rhinovirus",Ensemble)) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()
Rhinovirus_predictors <- read.csv(NA_S1_T2,sep=";") %>%
  dplyr::filter(grepl("Rhinovirus",Ensemble)) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(Rhinovirus_predictors)
Rhinovirus_predictors <- read.csv(NA_S2_T0,sep=";") %>%
  dplyr::filter(grepl("Rhinovirus",Ensemble)) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(Rhinovirus_predictors)
Rhinovirus_predictors <- read.csv(NA_S2_T2,sep=";") %>%
  dplyr::filter(grepl("Rhinovirus",Ensemble)) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(Rhinovirus_predictors)

## aydin
Rhinovirus_predictors <- read.csv(AY_S1_T0,sep=",") %>%
  dplyr::filter(grepl("Rhinovirus",VIRUS_TYPE)) %>%
  dplyr::select(PROBE_SET_ID) %>%
  dplyr::transmute(Predictor = PROBE_SET_ID) %>%
  probeIDs2GeneNames() %>%
  rbind(Rhinovirus_predictors)
Rhinovirus_predictors <- read.csv(AY_S1_T2,sep=",") %>%
  dplyr::filter(grepl("Rhinovirus",VIRUS_TYPE)) %>%
  dplyr::select(PROBE_SET_ID) %>%
  dplyr::transmute(Predictor = PROBE_SET_ID) %>%
  probeIDs2GeneNames() %>%
  rbind(Rhinovirus_predictors)
Rhinovirus_predictors <- read.csv(AY_S2_T0,sep=",") %>%
  dplyr::filter(grepl("Rhinovirus",VIRUS_TYPE)) %>%
  dplyr::select(PROBE_SET_ID) %>%
  dplyr::transmute(Predictor = PROBE_SET_ID) %>%
  probeIDs2GeneNames() %>%
  rbind(Rhinovirus_predictors)
Rhinovirus_predictors <- read.csv(AY_S2_T2,sep=",") %>%
  dplyr::filter(grepl("Rhinovirus",VIRUS_TYPE)) %>%
  dplyr::select(PROBE_SET_ID) %>%
  dplyr::transmute(Predictor = PROBE_SET_ID) %>%
  probeIDs2GeneNames() %>%
  rbind(Rhinovirus_predictors)
Rhinovirus_predictors <- read.csv(AY_S3_T0,sep=",",header = FALSE) %>%
  dplyr::filter(grepl("Rhinovirus",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(Rhinovirus_predictors)
Rhinovirus_predictors <- read.csv(AY_S3_T2,sep=",",header = FALSE) %>%
  dplyr::filter(grepl("Rhinovirus",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(Rhinovirus_predictors)

## SSN_Dream_Team
Rhinovirus_predictors <- read.csv(SS_S1_T0,sep=",",header=FALSE) %>%
  dplyr::filter(grepl("Rhinovirus",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(Rhinovirus_predictors)
Rhinovirus_predictors <- read.csv(SS_S1_T2,sep=",",header=FALSE) %>%
  dplyr::filter(grepl("Rhinovirus",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(Rhinovirus_predictors)
Rhinovirus_predictors <- read.csv(SS_S2_T0,sep=",",header=FALSE) %>%
  dplyr::filter(grepl("Rhinovirus",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(Rhinovirus_predictors)
Rhinovirus_predictors <- read.csv(SS_S2_T2,sep=",",header=FALSE) %>%
  dplyr::filter(grepl("Rhinovirus",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(Rhinovirus_predictors)
Rhinovirus_predictors <- read.csv(SS_S3_T0,sep=",",header=FALSE) %>%
  dplyr::filter(grepl("Rhinovirus",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(Rhinovirus_predictors)
Rhinovirus_predictors <- read.csv(SS_S3_T0,sep=",",header=FALSE) %>%
  dplyr::filter(grepl("Rhinovirus",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(Rhinovirus_predictors)

## Txsolo
Rhinovirus_predictors <- read.csv(TX_S1_T0,sep=",") %>%
  dplyr::filter(grepl("DUKE",virus)) %>%
  dplyr::select(feature) %>%
  dplyr::transmute(Predictor = feature) %>%
  unigeneIDs2GeneNames() %>%
  rbind(Rhinovirus_predictors)
Rhinovirus_predictors <- read.csv(TX_S1_T2,sep=",") %>%
  dplyr::filter(grepl("DUKE",virus)) %>%
  dplyr::select(feature) %>%
  dplyr::transmute(Predictor = feature) %>%
  unigeneIDs2GeneNames() %>%
  rbind(Rhinovirus_predictors)
Rhinovirus_predictors <- read.csv(TX_S2_T0,sep=",") %>%
  dplyr::filter(grepl("DUKE",virus)) %>%
  dplyr::select(feature) %>%
  dplyr::transmute(Predictor = feature) %>%
  unigeneIDs2GeneNames() %>%
  rbind(Rhinovirus_predictors)
Rhinovirus_predictors <- read.csv(TX_S2_T2,sep=",") %>%
  dplyr::filter(grepl("DUKE",virus)) %>%
  dplyr::select(feature) %>%
  dplyr::transmute(Predictor = feature) %>%
  unigeneIDs2GeneNames() %>%
  rbind(Rhinovirus_predictors)

## cwruPatho
Rhinovirus_predictors <- read.csv(CW_S1_T0,sep=",") %>%
  dplyr::filter(grepl("Rhinovirus",Virus)) %>%
  dplyr::select(Probe.Set.ID) %>%
  dplyr::transmute(Predictor = Probe.Set.ID) %>%
  probeIDs2GeneNames() %>%
  rbind(Rhinovirus_predictors)
Rhinovirus_predictors <- read.csv(CW_S1_T2,sep=",") %>%
  dplyr::filter(grepl("Rhinovirus",Virus)) %>%
  dplyr::select(Probe.Set.ID) %>%
  dplyr::transmute(Predictor = Probe.Set.ID) %>%
  probeIDs2GeneNames() %>%
  rbind(Rhinovirus_predictors)
Rhinovirus_predictors <- read.csv(CW_S2_T2,sep=",") %>%
  dplyr::filter(grepl("Rhinovirus",Virus)) %>%
  dplyr::select(Probe.Set.ID) %>%
  dplyr::transmute(Predictor = Probe.Set.ID) %>%
  probeIDs2GeneNames() %>%
  rbind(Rhinovirus_predictors)
Rhinovirus_predictors <- read.csv(CW_S3_T0,sep=",") %>%
  dplyr::filter(grepl("Rhinovirus",Virus)) %>%
  dplyr::select(Probe.Set.ID) %>%
  dplyr::transmute(Predictor = Probe.Set.ID) %>%
  probeIDs2GeneNames() %>%
  rbind(Rhinovirus_predictors)
Rhinovirus_predictors <- read.csv(CW_S3_T2,sep=",") %>%
  dplyr::filter(grepl("Rhinovirus",Virus)) %>%
  dplyr::select(Probe.Set.ID) %>%
  dplyr::transmute(Predictor = Probe.Set.ID) %>%
  probeIDs2GeneNames() %>%
  rbind(Rhinovirus_predictors)

## Aganita
Rhinovirus_predictors <- read.csv(AG_S1_T0,sep=",") %>%
  dplyr::filter(grepl("Rhinovirus",Study)) %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene) %>%
  rbind(Rhinovirus_predictors)
Rhinovirus_predictors <- read.csv(AG_S1_T2,sep=",") %>%
  dplyr::filter(grepl("Rhinovirus",Study)) %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene) %>%
  rbind(Rhinovirus_predictors)
Rhinovirus_predictors <- read.csv(AG_S2_T0,sep=",") %>%
  dplyr::filter(grepl("Rhinovirus",Study)) %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene) %>%
  rbind(Rhinovirus_predictors)
Rhinovirus_predictors <- read.csv(AG_S2_T2,sep=",") %>%
  dplyr::filter(grepl("Rhinovirus",Study)) %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene) %>%
  rbind(Rhinovirus_predictors)
Rhinovirus_predictors <- read.csv(AG_S3_T0,sep=",") %>%
  dplyr::filter(grepl("Rhinovirus",Study)) %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene) %>%
  rbind(Rhinovirus_predictors)
Rhinovirus_predictors <- read.csv(AG_S3_T2,sep=",") %>%
  dplyr::filter(grepl("Rhinovirus",Study)) %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene) %>%
  rbind(Rhinovirus_predictors)

# calculate support
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

plot(res, sort.by="p-value")

plot(res, Layout="landscape",degree=2:4, sort.by="p-value")

sigResDF <- summary(res)[["Table"]] %>%
  as.data.frame() %>%
  dplyr::filter(P.value < 0.005)

# end from https://cran.r-project.org/web/packages/SuperExactTest/vignettes/set_html.html

H1N1_predictors <- sigResDF %>%
  dplyr::select(Elements) %>%
  .[,] %>%
  strsplit(.,split=", ") %>%
  unlist() %>%
  unique()

H3N2_predictors <- H3N2_predictors %>%
  dplyr::filter(Predictor != "NA") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  table() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".",Support = "Freq") %>%
  dplyr::mutate(Predictor = as.character(Predictor),Support = as.numeric(Support))

Rhinovirus_predictors <- Rhinovirus_predictors %>%
  dplyr::filter(Predictor != "NA") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  table() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".",Support = "Freq") %>%
  dplyr::mutate(Predictor = as.character(Predictor),Support = as.numeric(Support))

H3N2_predictors %>%
  ggplot(aes(x=reorder(Support,-Support))) + 
  geom_bar(stat="count") +
  geom_text(stat="count",aes(label=..count..),vjust=-1,size=3) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylim(0,9000) +
  xlab("Gene Support (# Studies Including Gene)") +
  ylab("Gene Count (# Genes)") +
  ggtitle("H3N2 Predictor Support")

Rhinovirus_predictors %>%
  ggplot(aes(x=reorder(Support,-Support))) + 
  geom_bar(stat="count") +
  geom_text(stat="count",aes(label=..count..),vjust=-1,size=3) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylim(0,9000) +
  xlab("Gene Support (# Studies Including Gene)") +
  ylab("Gene Count (# Genes)") +
  ggtitle("Rhinovirus Predictor Support")

# reactome pathway enrichment analysis
H1N1_predictors %>%
  ReactomeBarplot(query_gene_symbols = .,
                  bar_color = "pink",
                  plot_title = "H1N1 Predictors With p-value < 0.005")

H3N2_predictors %>%
  dplyr::filter(Support >= 3) %>%
  dplyr::select(Predictor) %>%
  unlist() %>%
  as.character() %>%
  ReactomeBarplot(query_gene_symbols = .,
                  bar_color = "red",
                  plot_title = "H3N2 Predictors Supported by Two Models or More")

Rhinovirus_predictors %>%
  dplyr::filter(Support >= 3) %>%
  dplyr::select(Predictor) %>%
  unlist() %>%
  as.character() %>%
  ReactomeBarplot(query_gene_symbols = .,
                  bar_color = "gold",
                  plot_title = "Rhinovirus Predictors Supported by Two Models or More")
