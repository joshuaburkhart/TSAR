
library(magrittr)
library(dplyr)

table.s3 <- read.csv(file="~/tmp/Table_S3.csv",
                     header = TRUE, 
                     stringsAsFactors = FALSE) %>%
  dplyr::select(-Fisher.Exact.P,-Benjamini.Hochberg.Adj.P)
lit.table.full <- read.csv(file="~/tmp/LitGeneSetOverlapSignificance.csv",
                           header=TRUE,
                           stringsAsFactors = FALSE)

table.join <- table.s3 %>% dplyr::full_join(lit.table.full,
                                            by = c("Num.Neither",
                                                   "Num.RVD.Only",
                                                   "Num.LIT.Only",
                                                   "Num.Both",
                                                   "Overlapping.Genes"))

table.join <- table.join %>%
  dplyr::select(Comparison,RVD.List,Lit.Article,Lit.List,
                Num.Neither,Num.RVD.Only,Num.LIT.Only,Num.Both,
                Overlapping.Genes,Fisher.Exact.P,BH.Adj.P)

table.join %>% write.table(file="~/tmp/Table_S3.xlsx.csv",sep=",",row.names=FALSE)