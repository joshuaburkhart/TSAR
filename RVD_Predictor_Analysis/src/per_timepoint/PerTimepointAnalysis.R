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

JD_S3_T0 <- paste(DOWNLOAD_DIR, "Subchallenge3_jdn_Time0_Predictors.csv",sep="")
JD_S3_T2 <- paste(DOWNLOAD_DIR, "Subchallenge3_jdn_Time24_Predictors.csv",sep="")

JD_T0_predictors <- read.csv(JD_S3_T0,header = TRUE) %>%
  dplyr::transmute(Predictor = PREDICTOR) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()

JD_T2_predictors <- read.csv(JD_S3_T2,header = TRUE) %>%
  dplyr::transmute(Predictor = PREDICTOR) %>%
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

ER_S2_T0 <- paste(DOWNLOAD_DIR, "Subchallenge2_Espoir_Time0_Predictors.csv",sep="")
ER_S2_T2 <- paste(DOWNLOAD_DIR, "Subchallenge2_Espoir_Time24_Predictors.csv",sep="")
                 
ER_T0_predictors <- read.csv(ER_S2_T0,header=TRUE) %>%
  dplyr::transmute(Predictor = FEATUREID) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()

ER_T2_predictors <- read.csv(ER_S2_T2,header=FALSE) %>%
  dplyr::transmute(Predictor = V1) %>%
  dplyr::select(Predictor) %>%
  probeIDs2GeneNames()
 
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

# calculate support

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

CW_T0_predictors <- CW_T0_predictors[!duplicated(CW_T0_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

JD_T0_predictors <- JD_T0_predictors[!duplicated(JD_T0_predictors),] %>%
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

JH_T0_predictors <- JH_T0_predictors[!duplicated(JH_T0_predictors),] %>%
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

list.gene.sets <- c(
                    AY_T0_predictors,
                    SS_T0_predictors,
                    CW_T0_predictors,
                    ER_T0_predictors,
                    JD_T0_predictors,
                    JH_T0_predictors,
                    BK_T0_predictors)

names(list.gene.sets) <- c(
                           "AY",
                           "SS",
                           "CW",
                           "ER",
                           "JD",
                           "JH",
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

JD_T2_predictors <- JD_T2_predictors[!duplicated(JD_T2_predictors),] %>%
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

ER_T2_predictors <- ER_T2_predictors[!duplicated(ER_T2_predictors),] %>%
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

BK_T2_predictors <- BK_T2_predictors[!duplicated(BK_T2_predictors),] %>%
  as.character() %>%
  as.data.frame() %>%
  dplyr::transmute_(Predictor = ".") %>%
  dplyr::mutate(Predictor = as.character(Predictor)) %>%
  dplyr::filter(Predictor != "NA")

list.gene.sets <- c(
                    AY_T2_predictors,
                    SS_T2_predictors,
                    CW_T2_predictors,
                    ER_T2_predictors,
                    JD_T2_predictors,
                    JH_T2_predictors,
                    BK_T2_predictors)

names(list.gene.sets) <- c(
                           "AY",
                           "SS",
                           "CW",
                           "ER",
                           "JD",
                           "JH",
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
  write(paste(GENELIST_DIR,"T0_Predictors.txt",sep=""))
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
  write(paste(GENELIST_DIR,"T2_Predictors.txt",sep=""))
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
  write(paste(GENELIST_DIR,"T0_T2_Intersection.txt",sep=""))

T0_Only <- dplyr::setdiff(T0_predictors,T2_predictors)
T0_Only %>%
  write(paste(GENELIST_DIR,"T0_Only.txt",sep=""))

T2_Only <- dplyr::setdiff(T2_predictors,T0_predictors)
T2_Only %>%
  write(paste(GENELIST_DIR,"T2_Only.txt",sep=""))

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
