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
library(xlsxjars)
library(openxlsx)
library(ggplot2)
library(dplyr)
library(rJava)

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


# gloabal variables
P_VAL_THRESH = 0.005

ENTREZ_UNIVERSE <- AnnotationDbi::select(hgu133a2.db,keys=ls(hgu133a2ENTREZID),columns="ENTREZID") %>% .[,"ENTREZID"]
ENTREZ_UNIVERSE <- ENTREZ_UNIVERSE[!duplicated(ENTREZ_UNIVERSE)]

DOWNLOAD_DIR <- "../data/downloads/"
RESULTS_DIR <- "../results/"
FIGURES_DIR <- paste(RESULTS_DIR,"Figs/",sep="")
XLSX_DIR <- paste(RESULTS_DIR,"ExcelFiles/",sep="")

CT_S2_T0 <- paste(DOWNLOAD_DIR, "Subchallenge2_ChengzheTian_Time0_Predictors.csv",sep="")
CT_S2_T2 <- paste(DOWNLOAD_DIR, "Subchallenge2_ChengzheTian_Time24_Predictors.csv",sep="")
CT_S3_T0 <- paste(DOWNLOAD_DIR, "Subchallenge3_ChengzheTian_Time0_Predictors.csv",sep="")
CT_S3_T2 <- paste(DOWNLOAD_DIR, "Subchallenge3_ChengzheTian_Time24_Predictors.csv",sep="")
     
CT_T0_predictors <- read.csv(CT_S2_T0,header=FALSE,sep=",") %>%
  dplyr::transmute(Predictor = V1) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()
CT_T0_predictors <- read.csv(CT_S3_T0,header=FALSE,sep=",") %>%
  dplyr::transmute(Predictor = V1) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(CT_T0_predictors)

CT_T2_predictors <- read.csv(CT_S2_T2,header=FALSE,sep=",") %>%
  dplyr::transmute(Predictor = V1) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()
CT_T2_predictors <- read.csv(CT_S3_T2,header=FALSE,sep=",") %>%
  dplyr::transmute(Predictor = V1) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(CT_T2_predictors)
             
CF_S2_T0 <- paste(DOWNLOAD_DIR, "Subchallenge2_ChristoferF_Time0_Predictors.csv",sep="")
CF_S2_T2 <- paste(DOWNLOAD_DIR, "Subchallenge2_ChristoferF_Time24_Predictors.csv",sep="")
            
CF_T0_predictors <- read.csv(CF_S2_T0,header=FALSE,sep=",") %>%
  dplyr::transmute(Predictor = V1) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()

CF_T2_predictors <- read.csv(CF_S2_T2,header=FALSE,sep=",") %>%
  dplyr::transmute(Predictor = V1) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()
      
ES_S2_T0 <- paste(DOWNLOAD_DIR, "Subchallenge2_ES.SJ_PREDICTOMI_Time0_Predictors.csv",sep="")
ES_S2_T2 <- paste(DOWNLOAD_DIR, "Subchallenge2_ES.SJ_PREDICTOMI_Time24_Predictors.csv",sep="")
            
ES_T0_predictors <- read.csv(ES_S2_T0,header=TRUE,sep=",") %>%
  dplyr::transmute(Predictor = Variable) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()

ES_T2_predictors <- read.csv(ES_S2_T2,header=TRUE,sep=",") %>%
  dplyr::transmute(Predictor = Variable) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()
      
ER_S2_T0 <- paste(DOWNLOAD_DIR, "Subchallenge2_Espoir_Time0_Predictors.csv",sep="")
ER_S2_T2 <- paste(DOWNLOAD_DIR, "Subchallenge2_Espoir_Time24_Predictors.csv",sep="")
                 
ER_T0_predictors <- read.csv(ER_S2_T0,header=FALSE) %>%
  dplyr::transmute(Predictor = V1) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()

ER_T2_predictors <- read.csv(ER_S2_T2,header=FALSE) %>%
  dplyr::transmute(Predictor = V1) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()
 
GN_S2_T0 <- paste(DOWNLOAD_DIR, "Subchallenge2_Gustavssonlab_Nordlinglab_0_Predictors.csv",sep="")
GN_S2_T2 <- paste(DOWNLOAD_DIR, "Subchallenge2_Gustavssonlab_Nordlinglab_24_Predictors.csv",sep="")
GN_S3_T0 <- paste(DOWNLOAD_DIR, "Subchallenge3_Gustavssonlab_Nordlinglab_0_Predictors.csv",sep="")
GN_S3_T2 <- paste(DOWNLOAD_DIR, "Subchallenge3_Gustavssonlab_Nordlinglab_24_Predictors.csv",sep="")
           
GN_T0_predictors <- read.csv(GN_S2_T0,header=TRUE) %>%
  dplyr::transmute(Predictor = predictors_symptomatic_score_0h) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()
GN_T0_predictors <- read.csv(GN_S3_T0,header=TRUE) %>%
  dplyr::transmute(Predictor = predictors_logsymptomatic_0h) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(GN_T0_predictors)

GN_T2_predictors <- read.csv(GN_S2_T2,header=TRUE) %>%
  dplyr::transmute(Predictor = predictors_symptomatic_score_24h) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()
GN_T2_predictors <- read.csv(GN_S3_T2,header=TRUE) %>%
  dplyr::transmute(Predictor = predictors_logsymptomatic_24h) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(GN_T2_predictors)
       
JK_S2_T0 <- paste(DOWNLOAD_DIR, "Subchallenge2_JayHawks-RVDC_Time0_Predictors.csv",sep="")
JK_S2_T2 <- paste(DOWNLOAD_DIR, "Subchallenge2_JayHawks-RVDC_Time24_Predictors.csv",sep="")
JK_S3_T0 <- paste(DOWNLOAD_DIR, "Subchallenge3_JayHawks-RVDC_Time0_Predictors.csv",sep="")
JK_S3_T2 <- paste(DOWNLOAD_DIR, "Subchallenge3_JayHawks-RVDC_Time24_Predictors.csv",sep="")
            
JK_T0_predictors <- read.csv(JK_S2_T0,header=TRUE,sep=",") %>%
  dplyr::transmute(Predictor = PREDICTORS) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()
JK_T0_predictors <- read.csv(JK_S3_T0,header=TRUE,sep=",") %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(JK_T0_predictors)

JK_T2_predictors <- read.csv(JK_S2_T2,header=TRUE,sep=",") %>%
  dplyr::transmute(Predictor = PREDICTORS) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()
JK_T2_predictors <- read.csv(JK_S3_T2,header=TRUE,sep=",") %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(JK_T2_predictors)
      
RC_S2_T0 <- paste(DOWNLOAD_DIR, "Subchallenge2_Rchow_Time0_Predictors.csv",sep="")
RC_S2_T2 <- paste(DOWNLOAD_DIR, "Subchallenge2_Rchow_Time24_Predictors.csv",sep="")
RC_S3_T0 <- paste(DOWNLOAD_DIR, "Subchallenge3_Rchow_Time0_Predictors.csv",sep="")
RC_S3_T2 <- paste(DOWNLOAD_DIR, "Subchallenge3_Rchow_Time24_Predictors.csv",sep="")
            
RC_T0_predictors <- read.csv(RC_S2_T0,header=TRUE,sep=",") %>%
  dplyr::transmute(Predictor = Probe) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()
RC_T0_predictors <- read.csv(RC_S3_T0,header=TRUE,sep=",") %>%
  dplyr::transmute(Predictor = Probe) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(RC_T0_predictors)

RC_T2_predictors <- read.csv(RC_S2_T2,header=TRUE,sep=",") %>%
  dplyr::transmute(Predictor = Probe) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()
RC_T2_predictors <- read.csv(RC_S3_T2,header=TRUE,sep=",") %>%
  dplyr::transmute(Predictor = Probe) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(RC_T2_predictors)
      
SC_S2_T0 <- paste(DOWNLOAD_DIR, "Subchallenge2_SchrodingersCat_Time0_Predictors.csv",sep="")
SC_S2_T2 <- paste(DOWNLOAD_DIR, "Subchallenge2_SchrodingersCat_Time24_Predictors.csv",sep="")
SC_S3_T0 <- paste(DOWNLOAD_DIR, "Subchallenge3_SchrodingersCat_Time0_Predictors.csv",sep="")
SC_S3_T2 <- paste(DOWNLOAD_DIR, "Subchallenge3_SchrodingersCat_Time24_Predictors.csv",sep="")
                  
SC_T0_predictors <- read.csv(SC_S2_T0,header=TRUE) %>%
  dplyr::transmute(Predictor = x) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()
SC_T0_predictors <- read.csv(SC_S3_T0,header=TRUE) %>%
  dplyr::transmute(Predictor = x) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(SC_T0_predictors)

SC_T2_predictors <- read.csv(SC_S2_T2,header=TRUE) %>%
  dplyr::transmute(Predictor = x) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()
SC_T2_predictors <- read.csv(SC_S3_T2,header=TRUE) %>%
  dplyr::transmute(Predictor = x) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(SC_T2_predictors)

SH_S2_T0 <- paste(DOWNLOAD_DIR, "Subchallenge2_Shosty_UpToHour0_Predictors.csv",sep="")
SH_S2_T2 <- paste(DOWNLOAD_DIR, "Subchallenge2_Shosty_UpToHour24_Predictors.csv",sep="")
SH_S3_T0 <- paste(DOWNLOAD_DIR, "Subchallenge3_Shosty_UpToHour0_Predictors.csv",sep="")
SH_S3_T2 <- paste(DOWNLOAD_DIR, "Subchallenge3_Shosty_UpToHour24_Predictors.csv",sep="")
                 
SH_T0_predictors <- read.csv(SH_S2_T0,header=FALSE) %>%
  dplyr::transmute(Predictor = V1) %>%
  dplyr::select(Predictor)
SH_T0_predictors <- read.csv(SH_S3_T0,header=FALSE) %>%
  dplyr::transmute(Predictor = V1) %>%
  dplyr::select(Predictor) %>%
  rbind(SH_T0_predictors)

SH_T2_predictors <- read.csv(SH_S2_T2,header=FALSE) %>%
  dplyr::transmute(Predictor = V1) %>%
  dplyr::select(Predictor)
SH_T2_predictors <- read.csv(SH_S3_T2,header=FALSE) %>%
  dplyr::transmute(Predictor = V1) %>%
  dplyr::select(Predictor) %>%
  rbind(SH_T2_predictors)
 
SU_S2_T0 <- paste(DOWNLOAD_DIR, "Subchallenge2_SunilKumar_Time0_Predictors.csv",sep="")
SU_S2_T2 <- paste(DOWNLOAD_DIR, "Subchallenge2_SunilKumar_Time24_Predictors.csv",sep="")
SU_S3_T0 <- paste(DOWNLOAD_DIR, "Subchallenge3_SunilKumar_Time0_Predictors.csv",sep="")
SU_S3_T2 <- paste(DOWNLOAD_DIR, "Subchallenge3_SunilKumar_Time24_Predictors.csv",sep="")

SU_T0_predictors <- read.csv(SU_S2_T0,header=FALSE) %>%
  dplyr::transmute(Predictor = V1) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()
SU_T0_predictors <- read.csv(SU_S3_T0,header=FALSE) %>%
  dplyr::transmute(Predictor = V1) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(SU_T0_predictors)

SU_T2_predictors <- read.csv(SU_S2_T2,header=FALSE) %>%
  dplyr::transmute(Predictor = V1) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()
SU_T2_predictors <- read.csv(SU_S3_T2,header=FALSE) %>%
  dplyr::transmute(Predictor = V1) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(SU_T2_predictors)

TD_S2_T0 <- paste(DOWNLOAD_DIR, "Subchallenge2_TempleDABI_Time0_Predictors.csv",sep="")
TD_S2_T2 <- paste(DOWNLOAD_DIR, "Subchallenge2_TempleDABI_Time24_Predictors.csv",sep="")
TD_S3_T0 <- paste(DOWNLOAD_DIR, "Subchallenge3_TempleDABI_Time0_Predictors.csv",sep="")
TD_S3_T2 <- paste(DOWNLOAD_DIR, "Subchallenge3_TempleDABI_Time24_Predictors.csv",sep="")

TD_T0_predictors <- read.csv(TD_S2_T0,header=TRUE,sep=",") %>%
  dplyr::transmute(Predictor = Feature) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()
TD_T0_predictors <- read.csv(TD_S3_T0,header=TRUE,sep=",") %>%
  dplyr::transmute(Predictor = Feature) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(TD_T0_predictors)

TD_T2_predictors <- read.csv(TD_S2_T2,header=TRUE,sep=",") %>%
  dplyr::transmute(Predictor = Feature) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()
TD_T2_predictors <- read.csv(TD_S3_T2,header=TRUE,sep=",") %>%
  dplyr::transmute(Predictor = Feature) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(TD_T2_predictors)

VP_S2_T0 <- paste(DOWNLOAD_DIR, "Subchallenge2_ViResPred_Time0_Predictors.csv",sep="")
VP_S2_T2 <- paste(DOWNLOAD_DIR, "Subchallenge2_ViResPred_Time24_Predictors.csv",sep="")
VP_S3_T0 <- paste(DOWNLOAD_DIR, "Subchallenge3_ViResPred_Time0_Predictors.csv",sep="")
VP_S3_T2 <- paste(DOWNLOAD_DIR, "Subchallenge3_ViResPred_Time24_Predictors.csv",sep="")

VP_T0_predictors <- read.csv(VP_S2_T0,header = TRUE,sep=",") %>%
  dplyr::transmute(Predictor = predictors) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()
VP_T0_predictors <- read.csv(VP_S3_T0,header = FALSE,sep=",") %>%
  dplyr::transmute(Predictor = V2) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(VP_T0_predictors)

VP_T2_predictors <- read.csv(VP_S2_T2,header = TRUE,sep=",") %>%
  dplyr::transmute(Predictor = predictors) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()
VP_T2_predictors <- read.csv(VP_S3_T2,header = FALSE,sep=",") %>%
  dplyr::transmute(Predictor = V2) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(VP_T2_predictors)

BK_S3_T0 <- paste(DOWNLOAD_DIR, "Subchallenge3_burkhajo_Time0_Predictors.csv",sep="")
BK_S3_T2 <- paste(DOWNLOAD_DIR, "Subchallenge3_burkhajo_Time24_Predictors.csv",sep="")

BK_T0_predictors <- read.csv(BK_S3_T0, header = TRUE) %>%
  dplyr::transmute(Predictor = PROBE) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()
BK_T2_predictors <- read.csv(BK_S3_T2, header = TRUE) %>%
  dplyr::transmute(Predictor = PROBE) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()

HV_S2_T0 <- paste(DOWNLOAD_DIR, "Subchallenge2_hackvirus_Time0_Predictors.csv",sep="")
HV_S2_T2 <- paste(DOWNLOAD_DIR, "Subchallenge2_hackvirus_Time24_Predictors.csv",sep="")

HV_T0_predictors <- read.csv(HV_S2_T0,header=TRUE,sep="\t") %>%
  dplyr::transmute(Predictor = predictors) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()
HV_T2_predictors <- read.csv(HV_S2_T2,header=TRUE,sep="\t") %>%
  dplyr::transmute(Predictor = predictors) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()

JH_S2_T0 <- paste(DOWNLOAD_DIR, "Subchallenge2_jhou_Time0_Predictors.csv",sep="")
JH_S2_T2 <- paste(DOWNLOAD_DIR, "Subchallenge2_jhou_Time24_Predictors.csv",sep="")
JH_S3_T0 <- paste(DOWNLOAD_DIR, "Subchallenge3_jhou_Time0_Predictors.csv",sep="")
JH_S3_T2 <- paste(DOWNLOAD_DIR, "Subchallenge3_jhou_Time24_Predictors.csv",sep="")

JH_T0_predictors <- read.csv(JH_S2_T0,header = TRUE,sep=",") %>%
  dplyr::transmute(Predictor = X) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()
JH_T0_predictors <- read.csv(JH_S3_T0,header = TRUE,sep=",") %>%
  dplyr::transmute(Predictor = X) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(JH_T0_predictors)

JH_T2_predictors <- read.csv(JH_S2_T2,header = TRUE,sep=",") %>%
  dplyr::transmute(Predictor = X) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()
JH_T2_predictors <- read.csv(JH_S3_T2,header = TRUE,sep=",") %>%
  dplyr::transmute(Predictor = X) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames() %>%
  rbind(JH_T2_predictors)

NA_S2_T0 <- paste(DOWNLOAD_DIR,"Subchallenge2_Nautilus_Time0_Predictors.csv",sep="")
NA_S2_T2 <- paste(DOWNLOAD_DIR,"Subchallenge2_Nautilus_Time24_Predictors.csv",sep="")

NA_T0_predictors <- read.csv(NA_S2_T0,sep=";") %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()

NA_T2_predictors <- read.csv(NA_S2_T2,sep=";") %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()

AY_S2_T0 <- paste(DOWNLOAD_DIR,"Subchallenge2_aydin_Time0_Predictors.csv",sep="")
AY_S2_T2 <- paste(DOWNLOAD_DIR,"Subchallenge2_aydin_Time24_Predictors.csv",sep="")
AY_S3_T0 <- paste(DOWNLOAD_DIR,"Subchallenge3_aydin_Time0_Predictors.csv",sep="")
AY_S3_T2 <- paste(DOWNLOAD_DIR,"Subchallenge3_aydin_Time24_Predictors.csv",sep="")

AY_T0_predictors <- read.csv(AY_S2_T0,sep=",") %>%
  dplyr::select(PROBE_SET_ID) %>%
  dplyr::transmute(Predictor = PROBE_SET_ID) %>%
  probeIDs2GeneNames()
AY_T0_predictors <- read.csv(AY_S3_T0,sep=",",header = FALSE) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(AY_T0_predictors)

AY_T2_predictors <- read.csv(AY_S2_T2,sep=",") %>%
  dplyr::select(PROBE_SET_ID) %>%
  dplyr::transmute(Predictor = PROBE_SET_ID) %>%
  probeIDs2GeneNames()
AY_T2_predictors <- read.csv(AY_S3_T2,sep=",",header = FALSE) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(AY_T2_predictors)

SS_S2_T0 <- paste(DOWNLOAD_DIR,"Subchallenge2_SSN DREAM TEAM-0_Predictors.csv",sep="")
SS_S2_T2 <- paste(DOWNLOAD_DIR,"Subchallenge2_SSN DREAM TEAM-24_Predictors.csv",sep="")
SS_S3_T0 <- paste(DOWNLOAD_DIR,"Subchallenge3_SSN DREAM TEAM-0_Predictors.csv",sep="")
SS_S3_T2 <- paste(DOWNLOAD_DIR,"Subchallenge3_SSN DREAM TEAM-24_Predictors.csv",sep="")

SS_T0_predictors <- read.csv(SS_S2_T0,sep=",",header=FALSE) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames()
SS_T0_predictors <- read.csv(SS_S3_T0,sep=",",header=FALSE) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(SS_T0_predictors)

SS_T2_predictors <- read.csv(SS_S2_T2,sep=",",header=FALSE) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames()
SS_T2_predictors <- read.csv(SS_S3_T2,sep=",",header=FALSE) %>%
  dplyr::select(V2) %>%
  dplyr::transmute(Predictor = V2) %>%
  probeIDs2GeneNames() %>%
  rbind(SS_T2_predictors)

# Unigene IDs
# DEE4 = H1N1, DEE5 = H3N2, DUKE = Rhinovirus
TX_S2_T0 <- paste(DOWNLOAD_DIR,"SubChallenge2_TXsolo_Time0_Predictors.csv",sep="")
TX_S2_T2 <- paste(DOWNLOAD_DIR,"SubChallenge2_TXsolo_Time24_Predictors.csv",sep="")

TX_T0_predictors <- read.csv(TX_S2_T0,sep=",") %>%
  dplyr::select(feature) %>%
  dplyr::transmute(Predictor = feature) %>%
  unigeneIDs2GeneNames()

TX_T2_predictors <- read.csv(TX_S2_T2,sep=",") %>%
  dplyr::select(feature) %>%
  dplyr::transmute(Predictor = feature) %>%
  unigeneIDs2GeneNames()

CW_S2_T2 <- paste(DOWNLOAD_DIR,"Subchallenge2_cwruPatho_Time24_Predictors.csv",sep="")
CW_S3_T0 <- paste(DOWNLOAD_DIR,"Subchallenge3_cwruPatho_Time0_Predictors.csv",sep="")
CW_S3_T2 <- paste(DOWNLOAD_DIR,"Subchallenge3_cwruPatho_Time24_Predictors.csv",sep="")

CW_T0_predictors <- read.csv(CW_S3_T0,sep=",") %>%
  dplyr::select(Probe.Set.ID) %>%
  dplyr::transmute(Predictor = Probe.Set.ID) %>%
  probeIDs2GeneNames()

CW_T2_predictors <- read.csv(CW_S2_T2,sep=",") %>%
  dplyr::select(Probe.Set.ID) %>%
  dplyr::transmute(Predictor = Probe.Set.ID) %>%
  probeIDs2GeneNames()
CW_T2_predictors <- read.csv(CW_S3_T2,sep=",") %>%
  dplyr::select(Probe.Set.ID) %>%
  dplyr::transmute(Predictor = Probe.Set.ID) %>%
  probeIDs2GeneNames() %>%
  rbind(CW_T2_predictors)

#HGNC IDs
AG_S2_T0 <- paste(DOWNLOAD_DIR,"Subchallenge2_Aganita_Time0_Predictors.csv",sep="")
AG_S2_T2 <- paste(DOWNLOAD_DIR,"Subchallenge2_Aganita_Time24_Predictors.csv",sep="")
AG_S3_T0 <- paste(DOWNLOAD_DIR,"Subchallenge3_Aganita_Time0_Predictors.csv",sep="")
AG_S3_T2 <- paste(DOWNLOAD_DIR,"Subchallenge3_Aganita_Time24_Predictors.csv",sep="")

AG_T0_predictors <- read.csv(AG_S2_T0,sep=",") %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene)
AG_T0_predictors <- read.csv(AG_S3_T0,sep=",") %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene) %>%
  rbind(AG_T0_predictors)

AG_T2_predictors <- read.csv(AG_S2_T2,sep=",") %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene)
AG_T2_predictors <- read.csv(AG_S3_T2,sep=",") %>%
  dplyr::select(Gene) %>%
  dplyr::transmute(Predictor = Gene) %>%
  rbind(AG_T2_predictors)

# calculate support

NA_T0_predictors <- NA_T0_predictors[!duplicated(NA_T0_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

AY_T0_predictors <- AY_T0_predictors[!duplicated(AY_T0_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

SS_T0_predictors <- SS_T0_predictors[!duplicated(SS_T0_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

TX_T0_predictors <- TX_T0_predictors[!duplicated(TX_T0_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

CW_T0_predictors <- CW_T0_predictors[!duplicated(CW_T0_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

AG_T0_predictors <- AG_T0_predictors[!duplicated(AG_T0_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

CT_T0_predictors <- CT_T0_predictors[!duplicated(CT_T0_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

CF_T0_predictors <- CF_T0_predictors[!duplicated(CF_T0_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

ES_T0_predictors <- ES_T0_predictors[!duplicated(ES_T0_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

ER_T0_predictors <- ER_T0_predictors[!duplicated(ER_T0_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

GN_T0_predictors <- GN_T0_predictors[!duplicated(GN_T0_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

RC_T0_predictors <- RC_T0_predictors[!duplicated(RC_T0_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

SC_T0_predictors <- SC_T0_predictors[!duplicated(SC_T0_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

SH_T0_predictors <- SH_T0_predictors[!duplicated(SH_T0_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

SU_T0_predictors <- SU_T0_predictors[!duplicated(SU_T0_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

TD_T0_predictors <- TD_T0_predictors[!duplicated(TD_T0_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

VP_T0_predictors <- VP_T0_predictors[!duplicated(VP_T0_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

JH_T0_predictors <- JH_T0_predictors[!duplicated(JH_T0_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

JK_T0_predictors <- JK_T0_predictors[!duplicated(JK_T0_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

HV_T0_predictors <- HV_T0_predictors[!duplicated(HV_T0_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

BK_T0_predictors <- BK_T0_predictors[!duplicated(BK_T0_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

list.gene.sets <- c(#NA_T0_predictors,
                    AY_T0_predictors,
                    SS_T0_predictors,
                    #TX_T0_predictors,
                    CW_T0_predictors,
                    #AG_T0_predictors,
                    #CT_T0_predictors,
                    #CF_T0_predictors,
                    ES_T0_predictors,
                    #ER_T0_predictors,
                    #GN_T0_predictors,
                    #RC_T0_predictors,
                    #SC_T0_predictors,
                    #SH_T0_predictors,
                    #SU_T0_predictors,
                    #TD_T0_predictors,
                    #VP_T0_predictors,
                    #JH_T0_predictors,
                    #JK_T0_predictors,
                    #HV_T0_predictors,
                    BK_T0_predictors)

names(list.gene.sets) <- c(#"NA",
                           "AY",
                           "SS",
                           #"TX",
                           "CW",
                           #"AG",
                           #"CT",
                           #"CF",
                           "ES",
                           #"ER",
                           #"GN",
                           #"RC",
                           #"SC",
                           #"SH",
                           #"SU",
                           #"TD",
                           #"VP",
                           #"JH",
                           #"JK",
                           #"HV",
                           "BK")

# from https://cran.r-project.org/web/packages/SuperExactTest/vignettes/set_html.html
length.gene.sets <- sapply(list.gene.sets,
                           FUN = length)

total=ENTREZ_UNIVERSE %>% length() # 24515

num.expected.overlap=total*do.call(base::prod,as.list(length.gene.sets/total))

sapply(0:min(length.gene.sets),function(i) SuperExactTest::dpsets(i, length.gene.sets, n=total))

res=SuperExactTest::supertest(list.gene.sets, n=total)

svg(filename=paste(FIGURES_DIR,"T0_Model_Intersections_(Spiral).svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
plot(res, sort.by="p-value")
dev.off()

svg(filename=paste(FIGURES_DIR,"T0_Model_Intersections_(Bar_Graph).svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
plot(res, Layout="landscape",degree=2:4, sort.by="p-value")
dev.off()

sigResDF <- summary(res)[["Table"]] %>%
  as.data.frame() %>%
  dplyr::filter(P.value < P_VAL_THRESH)

# end from https://cran.r-project.org/web/packages/SuperExactTest/vignettes/set_html.html

T0_predictors <- sigResDF %>%
  dplyr::select(Elements) %>%
  .[,] %>%
  strsplit(.,split=", ") %>%
  unlist() %>%
  unique()

NA_T2_predictors <- NA_T2_predictors[!duplicated(NA_T2_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

AY_T2_predictors <- AY_T2_predictors[!duplicated(AY_T2_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

SS_T2_predictors <- SS_T2_predictors[!duplicated(SS_T2_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

TX_T2_predictors <- TX_T2_predictors[!duplicated(TX_T2_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

CW_T2_predictors <- CW_T2_predictors[!duplicated(CW_T2_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

AG_T2_predictors <- AG_T2_predictors[!duplicated(AG_T2_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

CT_T2_predictors <- CT_T2_predictors[!duplicated(CT_T2_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

CF_T2_predictors <- CF_T2_predictors[!duplicated(CF_T2_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

ES_T2_predictors <- ES_T2_predictors[!duplicated(ES_T2_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

ER_T2_predictors <- ER_T2_predictors[!duplicated(ER_T2_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

GN_T2_predictors <- GN_T2_predictors[!duplicated(GN_T2_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

RC_T2_predictors <- RC_T2_predictors[!duplicated(RC_T2_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

SC_T2_predictors <- SC_T2_predictors[!duplicated(SC_T2_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

SH_T2_predictors <- SH_T2_predictors[!duplicated(SH_T2_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

SU_T2_predictors <- SU_T2_predictors[!duplicated(SU_T2_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

TD_T2_predictors <- TD_T2_predictors[!duplicated(TD_T2_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

VP_T2_predictors <- VP_T2_predictors[!duplicated(VP_T2_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

JH_T2_predictors <- JH_T2_predictors[!duplicated(JH_T2_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

JK_T2_predictors <- JK_T2_predictors[!duplicated(JK_T2_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

HV_T2_predictors <- HV_T2_predictors[!duplicated(HV_T2_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

BK_T2_predictors <- BK_T2_predictors[!duplicated(BK_T2_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

list.gene.sets <- c(#NA_T2_predictors,
                    AY_T2_predictors,
                    SS_T2_predictors,
                    #TX_T2_predictors,
                    CW_T2_predictors,
                    #AG_T2_predictors,
                    #CT_T2_predictors,
                    #CF_T2_predictors,
                    ES_T2_predictors,
                    #ER_T2_predictors,
                    #GN_T2_predictors,
                    #RC_T2_predictors, #all predictors
                    #SC_T2_predictors, #all predictors
                    #SH_T2_predictors,
                    #SU_T2_predictors,
                    #TD_T2_predictors,
                    #VP_T2_predictors, #all predictors
                    #JH_T2_predictors,
                    #JK_T2_predictors,
                    #HV_T2_predictors,
                    BK_T2_predictors)

names(list.gene.sets) <- c(#"NA",
                           "AY",
                           "SS",
                           #"TX",
                           "CW",
                           #"AG",
                           #"CT",
                           #"CF",
                           "ES",
                           #"ER",
                           #"GN",
                           #"RC",
                           #"SC",
                           #"SH",
                           #"SU",
                           #"TD",
                           #"VP",
                           #"JH",
                           #"JK",
                           #"HV",
                           "BK")

# from https://cran.r-project.org/web/packages/SuperExactTest/vignettes/set_html.html
length.gene.sets <- sapply(list.gene.sets,
                           FUN = length)

total=ENTREZ_UNIVERSE %>% length() # 24515

num.expected.overlap=total*do.call(base::prod,as.list(length.gene.sets/total))

sapply(0:min(length.gene.sets),function(i) SuperExactTest::dpsets(i, length.gene.sets, n=total))

res=SuperExactTest::supertest(list.gene.sets, n=total)

svg(filename=paste(FIGURES_DIR,"T2_Model_Intersections_(Spiral).svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
plot(res, sort.by="p-value")
dev.off()

svg(filename=paste(FIGURES_DIR,"T2_Model_Intersections_(Bar_Graph).svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
plot(res, Layout="landscape",degree=2:4, sort.by="p-value")
dev.off()

sigResDF <- summary(res)[["Table"]] %>%
  as.data.frame() %>%
  dplyr::filter(P.value < P_VAL_THRESH)

# end from https://cran.r-project.org/web/packages/SuperExactTest/vignettes/set_html.html

T2_predictors <- sigResDF %>%
  dplyr::select(Elements) %>%
  .[,] %>%
  strsplit(.,split=", ") %>%
  unlist() %>%
  unique()

# reactome pathway enrichment analysis

## T0
title <- paste("T0 Predictors With p-value < ",P_VAL_THRESH,sep="")
T0_predictors %>%
  write.xlsx(paste(XLSX_DIR,"T0_Predictors.xlsx",sep=""))
svg(filename=paste(FIGURES_DIR,"T0_Predictor_PA.svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
T0_predictors %>%
  ReactomeBarplot(query_gene_symbols = .,
                  bar_color = "lightskyblue",
                  plot_title = title)
dev.off()

## T2
title <- paste("T2 Predictors With p-value < ",P_VAL_THRESH,sep="")
T2_predictors %>%
  write.xlsx(paste(XLSX_DIR,"T2_Predictors.xlsx",sep=""))
svg(filename=paste(FIGURES_DIR,"T2_Predictor_PA.svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
T2_predictors %>%
  ReactomeBarplot(query_gene_symbols = .,
                  bar_color = "midnightblue",
                  plot_title = title)
dev.off()

# Intersections

T0_T2_predictors <- SuperExactTest::intersect(T0_predictors,
                                                  T2_predictors)
T0_T2_predictors %>%
  write.xlsx(paste(XLSX_DIR,"T0_T2_Intersection.xlsx",sep=""))

T0_Only <- dplyr::setdiff(T0_predictors,T2_predictors)
T0_Only %>%
  write.xlsx(paste(XLSX_DIR,"T0_Only.xlsx",sep=""))

T2_Only <- dplyr::setdiff(T2_predictors,T0_predictors)
T2_Only %>%
  write.xlsx(paste(XLSX_DIR,"T2_Only.xlsx",sep=""))

svg(filename=paste(FIGURES_DIR,"PerTimepointVennDiagram.svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
VennDiagram::draw.pairwise.venn(area1=length(T0_predictors),
                              area2=length(T2_predictors),
                              category = c("T0","T2"),
                              cross.area=length(T0_T2_predictors),
                              fill=c("lightskyblue","midnightblue"))
dev.off()

## T0 T2
svg(filename=paste(FIGURES_DIR,"T0_T2_Predictor_PA.svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
T0_T2_predictors %>%
  ReactomeBarplot(query_gene_symbols = .,
                  bar_color = "blue",
                  plot_title = paste("T0 & T2 Predictors With p-value < ",P_VAL_THRESH))
dev.off()


## T0 only
svg(filename=paste(FIGURES_DIR,"T0-Only_Predictor_PA.svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
T0_Only %>%
  ReactomeBarplot(query_gene_symbols = .,
                  bar_color = "lightskyblue",
                  plot_title = paste("T0-Only Predictors with p-value < ",P_VAL_THRESH))
dev.off()

## T2 only
svg(filename=paste(FIGURES_DIR,"T2-Only_Predictor_PA.svg",sep=""),
    width=15,
    height=15,
    pointsize=12)
T2_Only %>%
  ReactomeBarplot(query_gene_symbols = .,
                  bar_color = "midnightblue",
                  plot_title = paste("T2-Only Predictors with p-value < ",P_VAL_THRESH))
dev.off()

system_command = paste("qlmanage -t -s 1000 -o ", FIGURES_DIR, " ",FIGURES_DIR,"*.svg",sep="")
system(paste(system_command))