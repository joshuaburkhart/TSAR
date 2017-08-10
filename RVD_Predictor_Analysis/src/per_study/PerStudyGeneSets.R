# ---
# title: Per Study Gene Sets
# author: Joshua Burkhart
# date: Aug 10, 2017
# ---

# gloabal variables
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

# load libraries
library(org.Hs.eg.db)
library(hgu133a2.db)
library(annotate)
library(magrittr)
library(dplyr)

# helper functions
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
  probeIDs <- unique(as.character(na.omit(annotate::lookUp(unigeneIDs, "hgu133a2.db", "UNIGENE"))))
  geneNames <- toupper(unique(as.character(na.omit(annotate::lookUp(probeIDs,"hgu133a2.db","SYMBOL")))))
  return(data.frame(Predictor = geneNames))
}


# H1N1

## Nautilus
H1N1_predictors <- read.csv(NA_S1_T0,sep=";") %>%
  dplyr::filter(grepl("H1N1",Ensemble)) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()
H1N1_predictors <- read.csv(NA_S1_T2,sep=";") %>%
  dplyr::filter(grepl("H1N1",Ensemble)) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(H1N1_predictors)
H1N1_predictors <- read.csv(NA_S2_T0,sep=";") %>%
  dplyr::filter(grepl("H1N1",Ensemble)) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(H1N1_predictors)
H1N1_predictors <- read.csv(NA_S2_T2,sep=";") %>%
  dplyr::filter(grepl("H1N1",Ensemble)) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(H1N1_predictors)

## aydin
H1N1_predictors <- read.csv(AY_S1_T0,sep=",") %>%
  dplyr::filter(grepl("H1N1",VIRUS_TYPE)) %>%
  dplyr::select(PROBE_SET_ID) %>%
  dplyr::transmute(Predictor = PROBE_SET_ID) %>%
  probeIDs2GeneNames() %>%
  rbind(H1N1_predictors)
H1N1_predictors <- read.csv(AY_S1_T2,sep=",") %>%
  dplyr::filter(grepl("H1N1",VIRUS_TYPE)) %>%
  dplyr::select(PROBE_SET_ID) %>%
  dplyr::transmute(Predictor = PROBE_SET_ID) %>%
  probeIDs2GeneNames() %>%
  rbind(H1N1_predictors)
H1N1_predictors <- read.csv(AY_S2_T0,sep=",") %>%
  dplyr::filter(grepl("H1N1",VIRUS_TYPE)) %>%
  dplyr::select(PROBE_SET_ID) %>%
  dplyr::transmute(Predictor = PROBE_SET_ID) %>%
  probeIDs2GeneNames() %>%
  rbind(H1N1_predictors)
H1N1_predictors <- read.csv(AY_S2_T2,sep=",") %>%
  dplyr::filter(grepl("H1N1",VIRUS_TYPE)) %>%
  dplyr::select(PROBE_SET_ID) %>%
  dplyr::transmute(Predictor = PROBE_SET_ID) %>%
  probeIDs2GeneNames() %>%
  rbind(H1N1_predictors)
H1N1_predictors <- read.csv(AY_S3_T0,sep=",",header = FALSE) %>%
  dplyr::filter(grepl("H1N1",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(H1N1_predictors)
H1N1_predictors <- read.csv(AY_S3_T2,sep=",",header = FALSE) %>%
  dplyr::filter(grepl("H1N1",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(H1N1_predictors)

## SSN_Dream_Team
H1N1_predictors <- read.csv(SS_S1_T0,sep=",",header=FALSE) %>%
  dplyr::filter(grepl("H1N1",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(H1N1_predictors)
H1N1_predictors <- read.csv(SS_S1_T2,sep=",",header=FALSE) %>%
  dplyr::filter(grepl("H1N1",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(H1N1_predictors)
H1N1_predictors <- read.csv(SS_S2_T0,sep=",",header=FALSE) %>%
  dplyr::filter(grepl("H1N1",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(H1N1_predictors)
H1N1_predictors <- read.csv(SS_S2_T2,sep=",",header=FALSE) %>%
  dplyr::filter(grepl("H1N1",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(H1N1_predictors)
H1N1_predictors <- read.csv(SS_S3_T0,sep=",",header=FALSE) %>%
  dplyr::filter(grepl("H1N1",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(H1N1_predictors)
H1N1_predictors <- read.csv(SS_S3_T0,sep=",",header=FALSE) %>%
  dplyr::filter(grepl("H1N1",V1)) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(H1N1_predictors)

## Txsolo
H1N1_predictors <- read.csv(TX_S1_T0,sep=",") %>%
  dplyr::filter(grepl("DEE4",virus)) %>%
  dplyr::select(feature) %>%
  dplyr::transmute(Predictor = feature) %>%
  unigeneIDs2GeneNames() %>%
  rbind(H1N1_predictors)
H1N1_predictors <- read.csv(TX_S1_T2,sep=",") %>%
  dplyr::filter(grepl("DEE4",virus)) %>%
  dplyr::select(feature) %>%
  dplyr::transmute(Predictor = feature) %>%
  unigeneIDs2GeneNames() %>%
  rbind(H1N1_predictors)
H1N1_predictors <- read.csv(TX_S2_T0,sep=",") %>%
  dplyr::filter(grepl("DEE4",virus)) %>%
  dplyr::select(feature) %>%
  dplyr::transmute(Predictor = feature) %>%
  unigeneIDs2GeneNames() %>%
  rbind(H1N1_predictors)
H1N1_predictors <- read.csv(TX_S2_T2,sep=",") %>%
  dplyr::filter(grepl("DEE4",virus)) %>%
  dplyr::select(feature) %>%
  dplyr::transmute(Predictor = feature) %>%
  unigeneIDs2GeneNames() %>%
  rbind(H1N1_predictors)

## cwruPatho
H1N1_predictors <- read.csv(CW_S1_T0,sep=",") %>%
  dplyr::filter(grepl("H1N1",Virus)) %>%
  dplyr::select(Probe.Set.ID) %>%
  dplyr::transmute(Predictor = Probe.Set.ID) %>%
  probeIDs2GeneNames() %>%
  rbind(H1N1_predictors)
H1N1_predictors <- read.csv(CW_S1_T2,sep=",") %>%
  dplyr::filter(grepl("H1N1",Virus)) %>%
  dplyr::select(Probe.Set.ID) %>%
  dplyr::transmute(Predictor = Probe.Set.ID) %>%
  probeIDs2GeneNames() %>%
  rbind(H1N1_predictors)
H1N1_predictors <- read.csv(CW_S2_T2,sep=",") %>%
  dplyr::filter(grepl("H1N1",Virus)) %>%
  dplyr::select(Probe.Set.ID) %>%
  dplyr::transmute(Predictor = Probe.Set.ID) %>%
  probeIDs2GeneNames() %>%
  rbind(H1N1_predictors)
H1N1_predictors <- read.csv(CW_S3_T0,sep=",") %>%
  dplyr::filter(grepl("H1N1",Virus)) %>%
  dplyr::select(Probe.Set.ID) %>%
  dplyr::transmute(Predictor = Probe.Set.ID) %>%
  probeIDs2GeneNames() %>%
  rbind(H1N1_predictors)
H1N1_predictors <- read.csv(CW_S3_T2,sep=",") %>%
  dplyr::filter(grepl("H1N1",Virus)) %>%
  dplyr::select(Probe.Set.ID) %>%
  dplyr::transmute(Predictor = Probe.Set.ID) %>%
  probeIDs2GeneNames() %>%
  rbind(H1N1_predictors)

## Aganita
H1N1_predictors <- read.csv(AG_S1_T0,sep=",") %>%
  dplyr::filter(grepl("H1N1",Study)) %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene) %>%
  rbind(H1N1_predictors)
H1N1_predictors <- read.csv(AG_S1_T2,sep=",") %>%
  dplyr::filter(grepl("H1N1",Study)) %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene) %>%
  rbind(H1N1_predictors)
H1N1_predictors <- read.csv(AG_S2_T0,sep=",") %>%
  dplyr::filter(grepl("H1N1",Study)) %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene) %>%
  rbind(H1N1_predictors)
H1N1_predictors <- read.csv(AG_S2_T2,sep=",") %>%
  dplyr::filter(grepl("H1N1",Study)) %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene) %>%
  rbind(H1N1_predictors)
H1N1_predictors <- read.csv(AG_S3_T0,sep=",") %>%
  dplyr::filter(grepl("H1N1",Study)) %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene) %>%
  rbind(H1N1_predictors)
H1N1_predictors <- read.csv(AG_S3_T2,sep=",") %>%
  dplyr::filter(grepl("H1N1",Study)) %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene) %>%
  rbind(H1N1_predictors)

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

# print sizes
H1N1_predictors %>% unique() %>% dim()
H3N2_predictors %>% unique() %>% dim()
Rhinovirus_predictors %>% unique() %>% dim()