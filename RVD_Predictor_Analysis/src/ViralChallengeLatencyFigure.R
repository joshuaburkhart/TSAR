
library(magrittr)
library(ggplot2)
library(dplyr)
#library(memisc)
#library(cowplot)

# Assume working directory matches source file location
# R Studio: Session->Set Working Directory->To Source File Location
#setwd("/media/burkhart/Media/Software/TSAR/src")
DATA_DIR <- "/media/burkhart/Media/Research/USB_TRANSFER/McWeeney/TSAR/"
TRAIN_SYMPTOM <- paste(DATA_DIR,"ViralChallenge_training_SymptomScoresByDay.tsv",sep="")
train_symptom_df <- read.table(TRAIN_SYMPTOM,header=T,sep="\t")
### Scan for Missing Data

sum(is.na(train_symptom_df)) #0
### Fix Missing Data

### Checking Distributions
#### Symptoms
# By Virus and Day

# RSV
train_symptom_df.s <- train_symptom_df %>% 
  dplyr::filter(STUDYID %>% as.character() == "DEE1 RSV")
train_symptom_df.m <- reshape2::melt(train_symptom_df.s[,c('STUDYDAY',
                                                           'SX_RUNNYNOSE',
                                                           'SX_COUGH',
                                                           'SX_HEADACHE',
                                                           'SX_MALAISE',
                                                           'SX_MYALGIA',
                                                           'SX_SNEEZE',
                                                           'SX_SORETHROAT',
                                                           'SX_STUFFYNOSE')],
                                     id.vars = 'STUDYDAY')

symptom_df.rsv <- data.frame(Virus = "RSV",
                             STUDYDAY = train_symptom_df.m$STUDYDAY,
                             variable = train_symptom_df.m$variable,
                             value = train_symptom_df.m$value)

# H1N1
train_symptom_df.s <- train_symptom_df %>% 
  dplyr::filter(STUDYID %>% as.character() %in%
                  c("DEE3 H1N1",
                    "DEE4X H1N1"))
train_symptom_df.m <- reshape2::melt(train_symptom_df.s[,c('STUDYDAY',
                                                           'SX_RUNNYNOSE',
                                                           'SX_COUGH',
                                                           'SX_HEADACHE',
                                                           'SX_MALAISE',
                                                           'SX_MYALGIA',
                                                           'SX_SNEEZE',
                                                           'SX_SORETHROAT',
                                                           'SX_STUFFYNOSE')],
                                     id.vars = 'STUDYDAY')

symptom_df.h1n1 <- data.frame(Virus = "H1N1",
                              STUDYDAY = train_symptom_df.m$STUDYDAY,
                              variable = train_symptom_df.m$variable,
                              value = train_symptom_df.m$value)

# H3N2
train_symptom_df.s <- train_symptom_df %>% 
  dplyr::filter(STUDYID %>% as.character() %in%
                  c("DEE2 H3N2",
                    "DEE5 H3N2"))
train_symptom_df.m <- reshape2::melt(train_symptom_df.s[,c('STUDYDAY',
                                                           'SX_RUNNYNOSE',
                                                           'SX_COUGH',
                                                           'SX_HEADACHE',
                                                           'SX_MALAISE',
                                                           'SX_MYALGIA',
                                                           'SX_SNEEZE',
                                                           'SX_SORETHROAT',
                                                           'SX_STUFFYNOSE')],
                                     id.vars = 'STUDYDAY')

symptom_df.h3n2 <- data.frame(Virus = "H3N2",
                              STUDYDAY = train_symptom_df.m$STUDYDAY,
                              variable = train_symptom_df.m$variable,
                              value = train_symptom_df.m$value)

# Rhinovirus
train_symptom_df.s <- train_symptom_df %>% 
  dplyr::filter(STUDYID %>% as.character() %in%
                  c("Rhinovirus Duke",
                    "Rhinovirus UVA"))
train_symptom_df.m <- reshape2::melt(train_symptom_df.s[,c('STUDYDAY',
                                                           'SX_RUNNYNOSE',
                                                           'SX_COUGH',
                                                           'SX_HEADACHE',
                                                           'SX_MALAISE',
                                                           'SX_MYALGIA',
                                                           'SX_SNEEZE',
                                                           'SX_SORETHROAT',
                                                           'SX_STUFFYNOSE')],
                                     id.vars = 'STUDYDAY')

symptom_df.rhinovirus <- data.frame(Virus = "Rhinovirus",
                                    STUDYDAY = train_symptom_df.m$STUDYDAY,
                                    variable = train_symptom_df.m$variable,
                                    value = train_symptom_df.m$value)

symptom_df <- rbind(symptom_df.rsv,
                    symptom_df.h1n1,
                    symptom_df.h3n2,
                    symptom_df.rhinovirus)

colnames(symptom_df) <- c("Virus","StudyDay","Symptom Name","SymptomValue")

symptom_df %>%
  dplyr::arrange(SymptomValue) %>%
  ggplot(aes(x=StudyDay,y=SymptomValue,fill=`Symptom Name`)) +
  geom_bar(stat="identity",position="stack") +
  facet_wrap( ~ Virus,ncol=1) +
  theme(strip.background = element_rect(fill="lightgray")) +
  coord_cartesian(ylim=c(0,125),xlim=c(-3,8)) +
  geom_vline(xintercept=0) +
  geom_vline(xintercept=1) +
  labs(x="Days Since Viral Challenge",
       y="Symptom Severity Scores") +
  scale_x_continuous(breaks=seq(-3,8,by=1))

