print("TSAR: Train 0h")
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
TRAIN_SYMPTOM <-
  paste(DATA_DIR,
        "ViralChallenge_training_SymptomScoresByDay.tsv",
        sep = "")
TRAIN_CLINICL <-
  paste(DATA_DIR, "ViralChallenge_training_CLINICAL.tsv", sep = "")
TRAIN_EXPRESS <-
  paste(DATA_DIR, "ViralChallenge_training_EXPRESSION_RMA.tsv", sep = "")

print("misc_config")
setwd(DATA_DIR)
cl <- makeCluster(detectCores()) 
registerDoParallel(cl)

print("helper functions")
rmse <- function(error)
{
  sqrt(mean(error^2))
}
add_virus_col <- function(x) {
  x %>%
    dplyr::mutate(VIRUS =
                    ifelse(
                      grepl(" H1N1$", STUDYID) %>% as.character(),
                      "H1N1",
                      ifelse(
                        grepl(" H3N2$", STUDYID) %>% as.character(),
                        "H3N2",
                        ifelse(
                          grepl(" RSV$", STUDYID) %>% as.character(),
                          "RSV",
                          ifelse(
                            grepl("^Rhinovirus ", STUDYID) %>% as.character(),
                            "Rhinovirus",
                            NA
                          )
                        )
                      )
                    )) %>%
    dplyr::mutate(VIRUS = VIRUS %>% as.factor()) %>%
    return()
}
add_tlt0_col <- function(x){
  x %>%
    dplyr::mutate(TLT0 = TIMEHOURS < 0) %>%
    return()
}

print("load_eset")
rma_exprs <- as.matrix(
  read.table(
    TRAIN_EXPRESS,
    header = TRUE,
    sep = "\t",
    row.names = 1,
    check.names = FALSE,
    as.is = TRUE
  )
)

print("load clinical data")
train_clinicl_df <- read.table(TRAIN_CLINICL, header = T, sep = "\t")

print("fix missing data")
train_clinicl_df <- train_clinicl_df %>%
  dplyr::mutate(EARLYTX = ifelse(is.na(EARLYTX), 0, 1))
train_clinicl_df <- train_clinicl_df %>%
  dplyr::mutate(SHAM = ifelse(is.na(SHAM), 0, 1))

print("keep only t<=0 arrays")
tlte0_cel <- train_clinicl_df %>%
  dplyr::filter(TIMEHOURS <= 0) %>%
  dplyr::mutate(CEL = as.character(CEL)) %>%
  dplyr::select(CEL)

tlte0_exprs <- rma_exprs[, tlte0_cel$CEL]

rma_eset <- ExpressionSet(assayData = tlte0_exprs)

print("filter_most_variant_genes")
top_2000 <- rma_eset %>%
  limma::lmFit() %>%
  limma::eBayes() %>%
  as.data.frame() %>%
  dplyr::add_rownames(var = "probesetid") %>%
  dplyr::top_n(2000, s2.post)

print("create_pheatmap")
if (IMAGES_AS_PNG) {
  png(
    filename = "/home/burkhart/Software/TSAR/src/rma_pheatmap.PNG",
    width = 20,
    height = 20,
    units="in",
    res=96,
    pointsize = 12
  )
}
rma_eset[top_2000$probesetid, ] %>% pheatmap()
if (IMAGES_AS_PNG) {
  dev.off()
}

print("snm")
cad.bio <- train_clinicl_df %>%
  dplyr::filter(TIMEHOURS <= 0) %>%
  add_tlt0_col() %>%
  dplyr::select(CEL,TLT0)
rownames(cad.bio) <- cad.bio$CEL
cad.bio <- cad.bio %>%
  dplyr::select(-CEL)

cad.adjrm = train_clinicl_df %>%
  dplyr::filter(TIMEHOURS <= 0) %>%
  add_tlt0_col() %>%
  dplyr::select(CEL,STUDYID,GENDER)
rownames(cad.adjrm) <- cad.adjrm$CEL
cad.adjrm <- cad.adjrm %>%
  dplyr::select(-CEL)

adj.var = model.matrix(~.,cad.adjrm)
bio.var = model.matrix(~.,cad.bio)
rma.data = rma_eset %>% exprs()

snmR.cad = snm(rma.data,
               bio.var,
               adj.var,
               rm.adj=TRUE,
               num.iter=1)

snm_eset <- snmR.cad$norm.dat %>% ExpressionSet()
  
print("filter_most_variant_genes")
top_2000 <- snm_eset %>%
  limma::lmFit() %>%
  limma::eBayes() %>%
  as.data.frame() %>%
  dplyr::add_rownames(var = "probesetid") %>%
  dplyr::top_n(2000, s2.post)

print("store features")
save(top_2000,file="/home/burkhart/Software/TSAR/src/top_2000.rda")

print("create_pheatmap")
if (IMAGES_AS_PNG) {
  png(
    filename = "/home/burkhart/Software/TSAR/src/snm_pheatmap.PNG",
    width = 20,
    height = 20,
    units="in",
    res=96,
    pointsize = 12
  )
}
snm_eset[top_2000$probesetid, ] %>% pheatmap()
if (IMAGES_AS_PNG) {
  dev.off()
}

print("SVM Regression LOOCV")

tlt0_exprs <- snm_eset[top_2000$probesetid, ] %>%
  exprs() %>%
  as.data.frame() %>%
  dplyr::select(match(train_clinicl_df %>%
                         dplyr::filter(TIMEHOURS < 0) %>%
                         dplyr::mutate(CEL = as.character(CEL)) %>%
                         dplyr::select(CEL) %>% .[,1],names(.)))

tlt0_scores <- train_clinicl_df[match(tlt0_exprs %>%
        colnames(),train_clinicl_df %>%
        dplyr::mutate(CEL = as.character(CEL)) %>%
        dplyr::select(CEL) %>% .[,1]),] %>%
  dplyr::mutate(SYMPTSCORE_SC3 = 10^LOGSYMPTSCORE_SC3) %>%
  dplyr::select(CEL, SYMPTSCORE_SC3)

training_data <- tlt0_exprs %>%
  t() %>%
  as.data.frame() %>%
  dplyr::add_rownames(var = "CEL") %>%
  dplyr::full_join(tlt0_scores,by="CEL") %>%
  dplyr::select(-CEL)

#trained_svm <- parallelSVM::parallelSVM(SYMPTSCORE_SC3 ~ .,
#                                  data = training_data,
#                                  numberCores=ifelse(detectCores() < 7,detectCores(),7),
#                                  cross = 10)

#tuneResult <- e1071::tune(e1071::svm, 
#                          SYMPTSCORE_SC3 ~ .,
#                          data = training_data,
#                          ranges = list(epsilon = seq(0,1,0.1), cost = 2^(2:9)),
#                          cross = 10)
#print(tuneResult)
#plot(tuneResult) # somewhere between 0.5 and 0.7

tuneResult <- e1071::tune(e1071::svm, 
                          SYMPTSCORE_SC3 ~ .,
                          data = training_data,
                          ranges = list(epsilon = seq(0.5,0.7,0.01), cost = 2^(2:9)),
                          cross = 10)
print(tuneResult)
plot(tuneResult)

print("test SVM")

tunedModel <- tuneResult$best.model
predicted <- predict(tunedModel,training_data) %>% as.numeric()
error <- training_data$SYMPTSCORE_SC3 - predicted
print(paste("RMSE:",rmse(error)))

print("store SVM")

save(trained_svm, file="/home/burkhart/Software/TSAR/src/trained_svm.rda")
