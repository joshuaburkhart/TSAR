
#libs
library(magrittr)
library(GeneOverlap)

#globs
DATA_DIR <- "/home/burkhart/Software/TSAR/RVD_Predictor_Analysis/data/genelists/"
RVD_DIR <- paste(DATA_DIR,"RVD/",sep="")
LIT_DIR <- paste(DATA_DIR,"LIT/",sep="")

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
    unique_lit_data[[paste(basename(filename2),column,sep="-")]] <- 
      unique(toupper(lit_data[[filename2]][,column]))
  }
}

#perform fisher exact test on columns
overlap_results <- list()
for(element1 in names(unique_rvd_data)){
  for(element2 in names(unique_lit_data)){
    overlap_results[[paste(element1,element2,sep="__")]] <- 
      newGeneOverlap(unique_rvd_data[[element1]],
    unique_lit_data[[element2]],
    spec="hg19.gene")
  }
}

#create table
genelist_comparison <- character()
fisher_exact_p <- numeric()
for(element in names(overlap_results)){
  genelist_comparison <- c(genelist_comparison,element)
  fisher_exact_p <- c(fisher_exact_p,getPval(testGeneOverlap(overlap_results[[element]])))
}
bh_adj_p <- p.adjust(fisher_exact_p,method="BH")
data.frame(Comparison = genelist_comparison,
                                     Fisher.Exact.P = fisher_exact_p,
                                     BH.Adj.P = bh_adj_p) %>%
  write.csv(file=paste(DATA_DIR,"LitGeneSetOverlapSignificance.csv"))
