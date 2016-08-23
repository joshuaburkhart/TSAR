# TSAR: Expression EDA
# author: Joshua Burkhart
# date: August 23, 2016

# libs
library(knitr)
library(magrittr)
library(dplyr)
library(affy)
library(ggplot2)
library(limma)
library(pheatmap)

# globals
IMAGES_AS_SVG <- TRUE
DATA_DIR <- "~/TSAR/"
TRAIN_SYMPTOM <- paste(DATA_DIR,"ViralChallenge_training_SymptomScoresByDay.tsv",sep="")
TRAIN_CLINICL <- paste(DATA_DIR,"ViralChallenge_training_CLINICAL.tsv",sep="")
PHAS1_CLINICL <- paste(DATA_DIR,"ViralChallenge_test_Phase1_CLINICAL.tsv",sep="")
PHAS2_CLINICL <- paste(DATA_DIR,"ViralChallenge_test_Phase2_CLINICAL.tsv",sep="")
PHAS3_CLINICL <- paste(DATA_DIR,"ViralChallenge_test_Phase3_CLINICAL.tsv",sep="")
TRAIN_EXPRESS <- paste(DATA_DIR,"ViralChallenge_training_EXPRESSION_RMA.tsv",sep="")
NS_ID_LIST <- c("3014",
                "3016",
                "3022",
                "3024",
                "5013",
                "flu003",
                "flu017",
                "13",
                "RSV004",
                "RSV019")

# misc_options
setwd(DATA_DIR) # for running from chunks
opts_knit$set(root.dir = DATA_DIR) # for knitting
par(mfrow=c(1,1)) # one row & one column per plot

# helper_functions
soi_label <- function(SHAM,LOGSYMPTSCORE_SC3,SUBJECTID,T25){
  ### sh: shams
  if(SHAM == 1){return("_sh")}
  ### jk: jackson scores > 60
  if(10^LOGSYMPTSCORE_SC3 > 60){return("_jk")}
  ### ns: no symptoms & no shedding
  if(SUBJECTID %in% NS_ID_LIST){return("_ns")}
  ### 25: top 25 quartile of jackson scores
  if(T25){return("_25")}
  return("")
}

# load_eset
rma_exprs <- as.matrix(read.table(TRAIN_EXPRESS,
                                  header=TRUE,
                                  sep="\t",
                                  row.names=1,
                                  as.is=TRUE))

# filter_CELs_of_interest
## load clinical data
train_clinicl_df <- read.table(TRAIN_CLINICL,header=T,sep="\t")

## fix missing data
train_clinicl_df <- train_clinicl_df %>%
  dplyr::mutate(EARLYTX = ifelse(is.na(EARLYTX),0,1))
train_clinicl_df <- train_clinicl_df %>% 
  dplyr::mutate(SHAM = ifelse(is.na(SHAM),0,1))

## calculate t<0 membership in top 25th percentile of jackson scores
top25 <- train_clinicl_df %>%
  dplyr::filter(TIMEHOURS < 0) %>%
  dplyr::filter(LOGSYMPTSCORE_SC3 > quantile(LOGSYMPTSCORE_SC3,.75)) %>%
  dplyr::select(SUBJECTID)

train_clinicl_df <- train_clinicl_df %>%
  dplyr::mutate(T25 = SUBJECTID %in% top25$SUBJECTID)

## label subjects of interest
train_clinicl_df <- train_clinicl_df %>%
  dplyr::rowwise() %>%
  dplyr::mutate(sid = paste(STUDYID,
                            "_",
                            SUBJECTID,
                            soi_label(SHAM,
                                      LOGSYMPTSCORE_SC3,
                                      SUBJECTID,
                                      T25),
                            sep=""))

## keep only t<0 arrays 
tlt0_sid_cel <- train_clinicl_df %>%
  dplyr::filter(TIMEHOURS < 0) %>%
  dplyr::select(sid,CEL)

tlt0_exprs <- exprs[,tlt0_sid_cel$CEL]

## rename matrix column names from CEL to sid
colnames(tlt0_exprs)[match(tlt0_sid_cel[,1],
        colnames(tlt0_exprs))] <- tlt0_sid_cel[,2][
    match(tlt0_sid_cel[,1],colnames(tlt0_exprs))]

rma_eset <- ExpressionSet(assayData=tlt0_exprs)

# filter_most_variant_genes
top_2000 <- rma_eset %>% 
  limma::lmFit(design=p_design) %>%
  limma::eBayes() %>%
  as.data.frame() %>%
  dplyr::add_rownames(var="probesetid") %>%
  dplyr::top_n(2000,s2.post)

# create_pheatmap
if(IMAGES_AS_SVG){
  svg(filename=paste(SRC_DIR,"affy_batch_pairwise_MAplot.svg",sep=""),
      width=15,
      height=15,
      pointsize=12)
}
rma_eset[top_2000$probesetid,] %>% pheatmap()
if(IMAGES_AS_SVG){
  dev.off()
}