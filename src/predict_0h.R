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
library(parallelSVM)
library(doParallel)

print("globals")
IMAGES_AS_PNG <- TRUE
DATA_DIR <- "/media/burkhart/D052-7853/McWeeney/TSAR/"
PHAS1_CLINICL <- paste(DATA_DIR,"ViralChallenge_test_Phase1_CLINICAL.tsv",sep="")
PHAS1_EXPRESS <- paste(DATA_DIR,"ViralChallenge_test_Phase1_EXPRESSION_RMA.tsv",sep="")

print("misc_config")
setwd(DATA_DIR)
cl <- makeCluster(detectCores()) 
registerDoParallel(cl)

# load features
load("/home/burkhart/Software/TSAR/src/top_2000.rda") #loads top_2000 from train_0h

# load trained model
load("/home/burkhart/Software/TSAR/src/trained_svm.rda") #loads trained_svm from train_0h

