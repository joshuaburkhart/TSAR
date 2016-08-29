print("TSAR: Predict 0h")
print(" author: Joshua Burkhart")
print(" date: August 27, 2016")

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
PHAS1_CLINICL <- paste(DATA_DIR,"ViralChallenge_test_Phase1_CLINICAL.tsv",sep="")
PHAS1_EXPRESS <- paste(DATA_DIR,"ViralChallenge_test_Phase1_EXPRESSION_RMA.tsv",sep="")

print("misc_config")
setwd(DATA_DIR)
cl <- makeCluster(detectCores()) 
registerDoParallel(cl)

print("load features")
load("/home/burkhart/Software/TSAR/src/top_2000.rda") #loads top_2000 from train_0h

print("load trained model")
load("/home/burkhart/Software/TSAR/src/tunedModel.rda") #loads trained_svm from train_0h

print("load phas1 clinical data")
phas1_clinicl_df <- read.table(PHAS1_CLINICL, header = T, sep = "\t")

print("load phas1 eset")
rma_exprs <- as.matrix(
  read.table(
    PHAS1_EXPRESS,
    header = TRUE,
    sep = "\t",
    row.names = 1,
    check.names = FALSE,
    as.is = TRUE
  )
)

print("keep only t<0 arrays")
tlt0_cel <- phas1_clinicl_df %>%
  dplyr::filter(TIMEHOURS < 0) %>%
  dplyr::mutate(CEL = as.character(CEL)) %>%
  dplyr::select(CEL)

tlt0_exprs <- rma_exprs[, tlt0_cel$CEL]

rma_eset <- ExpressionSet(assayData = tlt0_exprs)

print("keep only training features")
tlt0_exprs <- rma_eset[top_2000$probesetid, ] %>%
  exprs() %>%
  as.data.frame() %>%
  dplyr::select(match(phas1_clinicl_df %>%
                        dplyr::filter(TIMEHOURS < 0) %>%
                        dplyr::mutate(CEL = as.character(CEL)) %>%
                        dplyr::select(CEL) %>% .[,1],names(.)))

phas1_data <- tlt0_exprs %>%
  t() %>%
  as.data.frame() %>%
  dplyr::add_rownames(var = "CEL") %>%
  dplyr::select(-CEL)

print("predict phas1 scores")
phas1_predictions <- predict(tunedModel,phas1_data) %>% as.numeric()

print("store phas1 predictions")
write.table(phas1_predictions,
            file="/home/burkhart/Software/TSAR/src/phas1_predictions.csv",
            sep=",")


