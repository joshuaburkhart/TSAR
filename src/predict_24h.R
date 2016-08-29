print("TSAR: Predict 24h")
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
DATA_DIR <- "/media/burkhart/D052-7853/McWeeney/TSAR/"
PHAS3_CLINICL <- paste(DATA_DIR,"ViralChallenge_test_Phase3_CLINICAL.tsv",sep="")
PHAS3_EXPRESS <- paste(DATA_DIR,"ViralChallenge_test_Phase3_EXPRESSION_RMA.tsv",sep="")

print("misc_config")
setwd(DATA_DIR)

print("helper functions")
filter_24h_times <- function(x){
  x %>% dplyr::filter(TIMEHOURS == 0 | # all studies
                        TIMEHOURS == 21.5 | # DEE1, 2, 3, 4X, 
                        TIMEHOURS == 18 | # DEE5
                        TIMEHOURS == 24) %>% # Rhinovirus Duke, UVA
    return()
}

print("load features")
load("/home/burkhart/Software/TSAR/src/top_2000.rda") #loads top_2000 from train_24h

print("load trained model")
load("/home/burkhart/Software/TSAR/src/tunedModel.rda") #loads trained_svm from train_24h

print("load phas3 clinical data")
phas3_clinicl_df <- read.table(PHAS3_CLINICL, header = T, sep = "\t")

print("load phas3 eset")
rma_exprs <- as.matrix(
  read.table(
    PHAS3_EXPRESS,
    header = TRUE,
    sep = "\t",
    row.names = 1,
    check.names = FALSE,
    as.is = TRUE
  )
)

print("keep only t<0 arrays")
ta24_cel <- phas3_clinicl_df %>%
  filter_24h_times() %>%
  dplyr::mutate(CEL = as.character(CEL)) %>%
  dplyr::select(CEL)

ta24_exprs <- rma_exprs[, ta24_cel$CEL]

rma_eset <- ExpressionSet(assayData = ta24_exprs)

print("keep only training features")
ta24_exprs <- rma_eset[top_2000$probesetid, ] %>%
  exprs() %>%
  as.data.frame() %>%
  dplyr::select(match(phas3_clinicl_df %>%
                        filter_24h_times() %>%
                        dplyr::mutate(CEL = as.character(CEL)) %>%
                        dplyr::select(CEL) %>% .[,1],names(.)))

phas3_data <- ta24_exprs %>%
  t() %>%
  as.data.frame() %>%
  dplyr::add_rownames(var = "CEL") %>%
  dplyr::select(-CEL)

print("predict phas3 scores")
phas3_predictions <- predict(tunedModel,phas3_data) %>% as.numeric()

print("store phas3 predictions")
write.table(phas3_predictions,
            file="/home/burkhart/Software/TSAR/src/phas3_predictions.csv",
            sep=",")


