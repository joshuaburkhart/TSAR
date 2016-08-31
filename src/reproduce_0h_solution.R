set.seed(88)
start_time <- proc.time()

print("Reproduce 0h Solution")
print("author: Joshua Burkhart")
print("date: August 30, 2016")
print("start time:")
print(start_time)

OUT_DIR <- "/home/burkhart/Software/TSAR/src/"
DATA_DIR <- "/media/burkhart/D052-7853/McWeeney/TSAR/"
TRAIN_CLINICL <- paste(DATA_DIR, "ViralChallenge_training_CLINICAL.tsv", sep = "")
TRAIN_EXPRESS <- paste(DATA_DIR, "ViralChallenge_training_EXPRESSION_RMA.tsv", sep = "")
PHAS1_CLINICL <- paste(DATA_DIR,"ViralChallenge_test_Phase1_CLINICAL.tsv",sep="")
PHAS1_EXPRESS <- paste(DATA_DIR,"ViralChallenge_test_Phase1_EXPRESSION_RMA.tsv",sep="")
PHAS2_CLINICL <- paste(DATA_DIR,"ViralChallenge_test_Phase2_CLINICAL.tsv",sep="")
PHAS2_EXPRESS <- paste(DATA_DIR,"ViralChallenge_test_Phase2_EXPRESSION_RMA.tsv",sep="")
PHAS3_CLINICL <- paste(DATA_DIR,"ViralChallenge_test_Phase3_CLINICAL.tsv",sep="")
PHAS3_EXPRESS <- paste(DATA_DIR,"ViralChallenge_test_Phase3_EXPRESSION_RMA.tsv",sep="")

setwd(DATA_DIR)

# load training rma expression data
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

# load training clinical data
train_clinicl_df <- read.table(TRAIN_CLINICL, header = T, sep = "\t")

# remove sham subjects and data after 24h
train_clinicl_df <- train_clinicl_df %>%
  dplyr::mutate(EARLYTX = ifelse(is.na(EARLYTX), 0, 1)) %>%
  dplyr::filter(is.na(SHAM)) %>%
  dplyr::filter(TIMEHOURS <= 24)

# select each subject's earliest and latest microarrays
train_earliest <- train_clinicl_df %>%
  dplyr::group_by(SUBJECTID) %>%
  dplyr::arrange(TIMEHOURS) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(MICROARRAY_TIME = as.factor("early"))

train_latest <- train_clinicl_df %>%
  dplyr::group_by(SUBJECTID) %>%
  dplyr::arrange(TIMEHOURS) %>%
  dplyr::slice(n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(MICROARRAY_TIME = as.factor("late"))

train_earliest_latest <- rbind(train_earliest,train_latest) %>% as.data.frame()

train_earliest_latest_exprs <- train_earliest_latest %>%
  dplyr::mutate(CEL = as.character(CEL)) %>%
  dplyr::select(CEL)

train_earliest_latest_eset <- ExpressionSet(
  assayData = rma_exprs[, train_earliest_latest_exprs$CEL])

# snm
cad.bio <- train_earliest_latest %>%
  dplyr::select(CEL,MICROARRAY_TIME)
rownames(cad.bio) <- cad.bio$CEL
cad.bio <- cad.bio %>%
  dplyr::select(-CEL)

cad.adjrm = train_earliest_latest %>%
  dplyr::select(CEL,STUDYID,GENDER)
rownames(cad.adjrm) <- cad.adjrm$CEL
cad.adjrm <- cad.adjrm %>%
  dplyr::select(-CEL)

adj.var = model.matrix(~.,cad.adjrm)
bio.var = model.matrix(~.,cad.bio)
rma.data = train_earliest_latest_eset %>% exprs()

snmR.cad = snm(rma.data,
               bio.var,
               adj.var,
               rm.adj=TRUE,
               num.iter=5)

snm_eset <- ExpressionSet(assayData = snmR.cad$norm.dat)

# select most variable probes
top_2000 <- snm_eset %>%
  limma::lmFit() %>%
  limma::eBayes() %>%
  as.data.frame() %>%
  dplyr::add_rownames(var = "probesetid") %>%
  dplyr::filter(!grepl("AFFX",probesetid)) %>%
  dplyr::top_n(2000, s2.post)

print("List of predictors used in model")
predictor_list <- data.frame(
  PROBE=top_2000$probesetid,
  COEF=top_2000$coefficients
) %>%
  dplyr::arrange(desc(COEF)) %>%
  dplyr::select(PROBE)

write.table(predictor_list,
            file=paste(OUT_DIR,"Subchallenge3_burkhajo_Time0_Predictors.csv",sep=""),
            sep=",",
            row.names=FALSE,
            quote=FALSE)

# select each subject's earliest microarray and power log score for svm training
earliest_exprs <- snm_eset[top_2000$probesetid, ] %>%
  exprs() %>%
  as.data.frame() %>%
  dplyr::select(match(train_earliest_latest %>%
                        dplyr::filter(as.character(MICROARRAY_TIME) == "early") %>%
                        dplyr::mutate(CEL = as.character(CEL)) %>%
                        dplyr::select(CEL) %>% .[,1],names(.)))

earliest_scores <- train_earliest_latest[match(earliest_exprs %>%
                                        colnames(),train_earliest_latest %>%
                                        dplyr::mutate(CEL = as.character(CEL)) %>%
                                        dplyr::select(CEL) %>% .[,1]),] %>%
  dplyr::mutate(SYMPTSCORE_SC3 = 10^LOGSYMPTSCORE_SC3) %>%
  dplyr::select(CEL, SUBJECTID, SYMPTSCORE_SC3)

earliest_training_data <- earliest_exprs %>%
  t() %>%
  as.data.frame() %>%
  dplyr::add_rownames(var = "CEL") %>%
  dplyr::full_join(earliest_scores,by="CEL")

# tune svm using k=10-fold cross validation
tunedModels <- e1071::tune(e1071::svm, 
                          SYMPTSCORE_SC3 ~ .,
                          data = earliest_training_data %>%
                            dplyr::select(-CEL, -SUBJECTID) ,
                          ranges = list(epsilon = seq(0,1,0.1), cost = 2^(2:9)),
                          cross = 10)
best_svm <- tunedModels$best.model

print("Leave-one-out Cross-Validations in the Training Data")
earliest_loocv <- data.frame(
  SUBJECTID=earliest_training_data$SUBJECTID,
  LOGSYMPTSCORE_SC3=NA)

for(out_row in 1:nrow(earliest_training_data)){
  print(paste("holding out row ",out_row,"...",sep=""))
  
  loo_earliest_test_row <- earliest_training_data %>%
    dplyr::select(-CEL,-SUBJECTID) %>%
    .[out_row,]
  
  earliest_loocv[out_row,2] <- predict(best_svm, loo_earliest_test_row) %>% log10()
}

write.table(earliest_loocv,
            file=paste(OUT_DIR,"Subchallenge3_burkhajo_Time0_LOOCVs.csv",sep=""),
            sep=",",
            row.names=FALSE,
            quote=FALSE)

print("Leaderboard Predictions")
# load testing rma expression data
phas1_rma_exprs <- as.matrix(
  read.table(
    PHAS1_EXPRESS,
    header = TRUE,
    sep = "\t",
    row.names = 1,
    check.names = FALSE,
    as.is = TRUE
  )
)

phas2_rma_exprs <- as.matrix(
  read.table(
    PHAS2_EXPRESS,
    header = TRUE,
    sep = "\t",
    row.names = 1,
    check.names = FALSE,
    as.is = TRUE
  )
)

phas3_rma_exprs <- as.matrix(
  read.table(
    PHAS3_EXPRESS,
    header = TRUE,
    sep = "\t",
    row.names = 1,
    check.names = FALSE,
    as.is = TRUE
  )
)

test_rma_eset <- ExpressionSet(
  assayData = rbind(phas1_rma_exprs %>% t(),
                        phas2_rma_exprs %>% t(),
                        phas3_rma_exprs %>% t()) %>% t())

# load testing clinical data
phas1_clinicl_df <- read.table(PHAS1_CLINICL, header = T, sep = "\t")
phas2_clinicl_df <- read.table(PHAS2_CLINICL, header = T, sep = "\t")
phas3_clinicl_df <- read.table(PHAS3_CLINICL, header = T, sep = "\t")

test_clinicl_df <- rbind(phas1_clinicl_df,
                         phas2_clinicl_df,
                         phas3_clinicl_df)

# remove sham subjects and data after 24h
test_clinicl_df <- test_clinicl_df %>%
  dplyr::mutate(EARLYTX = ifelse(is.na(EARLYTX), 0, 1)) %>%
  dplyr::filter(is.na(SHAM)) %>%
  dplyr::filter(TIMEHOURS <= 24)

# select each subject's earliest microarray
test_earliest <- test_clinicl_df %>%
  dplyr::group_by(SUBJECTID) %>%
  dplyr::arrange(TIMEHOURS) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  as.data.frame()

test_earliest_exprs <- test_rma_eset[top_2000$probesetid, ] %>%
  exprs() %>%
  as.data.frame() %>%
  dplyr::select(match(test_earliest %>%
                        dplyr::mutate(CEL = as.character(CEL)) %>%
                        dplyr::select(CEL) %>% .[,1],names(.)))

test_earliest_subjects <- test_earliest[match(test_earliest_exprs %>%
                                                 colnames(),test_earliest %>%
                                                 dplyr::mutate(CEL = as.character(CEL)) %>%
                                                 dplyr::select(CEL) %>% .[,1]),] %>%
                                              dplyr::select(CEL, SUBJECTID)

test_earliest_data <- test_earliest_exprs %>%
  t() %>%
  as.data.frame() %>%
  dplyr::add_rownames(var = "CEL") %>%
  dplyr::full_join(test_earliest_subjects,by="CEL")

# predict scores
test_predictions <- predict(best_svm,
                            test_earliest_data %>%
                              dplyr::select(-CEL,-SUBJECTID)) %>%
                             log10()

leaderboard_submission <- data.frame(
  SUBJECTID = test_earliest_data$SUBJECTID,
  LOGSYMPTSCORE_SC3 = test_predictions) %>%
  dplyr::arrange(SUBJECTID)

write.table(leaderboard_submission,
            file="/home/burkhart/Software/TSAR/src/Subchallenge3_burkhajo_Time0_LeaderboardPredictions.csv",
            sep=",",
            row.names=FALSE,
            quote=FALSE)

print("end time:")
print(proc.time() - start_time)
