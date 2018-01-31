
#libs
library(magrittr)
library(hgu133a2.db)
library(GeneOverlap)

#globs
DATA_DIR <- "/home/burkhart/Software/TSAR/RVD_Predictor_Analysis/data/genelists/"
RVD_DIR <- paste(DATA_DIR,"RVD/",sep="")
LIT_DIR <- paste(DATA_DIR,"LIT/",sep="")
OVERLAP_DATA_TABLE_FN <- paste(DATA_DIR,"LitGeneSetOverlapSignificance.csv",sep="")
OVERLAP_LISTS_FN <- paste(DATA_DIR,"LitGeneSetIntersection.csv",sep="")
GENE_UNIVERSE_SIZE <- as.list(hgu133a2SYMBOL[hgu133a2SYMBOL %>% mappedkeys()]) %>% unlist(use.names = FALSE) %>% unique() %>% length()

#load RVD Challenge gene lists (single columns without headers)
rvd_filenames <- list.files(RVD_DIR, pattern="*.csv", full.names=TRUE)
rvd_data <- lapply(rvd_filenames, read.csv, header=FALSE, stringsAsFactors=FALSE)
names(rvd_data) <- rvd_filenames

#load gene lists from literature (1+ columns with headers)
lit_filenames <- list.files(LIT_DIR, pattern="*.csv", full.names=TRUE)
lit_data <- lapply(lit_filenames, read.csv, header=TRUE, stringsAsFactors=FALSE)
names(lit_data) <- lit_filenames

#store unique column values
unique_rvd_data <- list()
for(filename1 in rvd_filenames){
  unique_rvd_data[[basename(filename1)]] <- unique(toupper(rvd_data[[filename1]][,1]))
}
unique_lit_data <- list()
for(filename2 in lit_filenames){
  for(column in colnames(lit_data[[filename2]])){
    unique_lit_data[[paste(basename(filename2),column,sep=",")]] <- 
      unique(toupper(lit_data[[filename2]][,column]))
  }
}

#perform fisher exact test on columns
overlap_results <- list()
for(element1 in names(unique_rvd_data)){
  for(element2 in names(unique_lit_data)){
    overlap_results[[paste(element1,element2,sep=",")]] <- 
      newGeneOverlap(
        unique_rvd_data[[element1]], # List A
        unique_lit_data[[element2]], # List B
        genome.size=GENE_UNIVERSE_SIZE)
  }
}

#create table
genelist_comparison <- character()
fisher_exact_p <- numeric()
contingency_table_notB_notA <- numeric()
contingency_table_notB_inA <- numeric()
contingency_table_inB_notA <- numeric()
contingency_table_inB_inA <- numeric()
intersection_list <- list()
for(element in names(overlap_results)){
  gene_overlap <- testGeneOverlap(overlap_results[[element]])
  cont_table <- getContbl(gene_overlap)
  
  genelist_comparison <- c(genelist_comparison,element)
  fisher_exact_p <- c(fisher_exact_p,getPval(gene_overlap))
  contingency_table_notB_notA <- c(contingency_table_notB_notA, cont_table['notB','notA'])
  contingency_table_notB_inA <- c(contingency_table_notB_inA, cont_table['notB','inA'])
  contingency_table_inB_notA <- c(contingency_table_inB_notA, cont_table['inB','notA'])
  contingency_table_inB_inA <- c(contingency_table_inB_inA, cont_table['inB','inA'])
  
  intersection_list[[element]] <- c(element,getIntersection(gene_overlap))
}
bh_adj_p <- p.adjust(fisher_exact_p,method="BH")
data.frame(Comparison = genelist_comparison,
           Num.Neither = contingency_table_notB_notA,
           Num.RVD.Only = contingency_table_notB_inA,
           Num.LIT.Only = contingency_table_inB_notA,
           Num.Both = contingency_table_inB_inA,
           Fisher.Exact.P = fisher_exact_p,
           BH.Adj.P = bh_adj_p) %>%
  write.csv(file=OVERLAP_DATA_TABLE_FN)

writeLines(unlist(lapply(intersection_list, paste, collapse=" ")),con = OVERLAP_LISTS_FN)


