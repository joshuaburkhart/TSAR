library(magrittr)
library(dplyr)

data.dir <- "/home/burkhart/Software/TSAR/RVD_Predictor_Analysis/data/"
genelist.dir <- "/home/burkhart/Software/TSAR/RVD_Predictor_Analysis/results/GeneLists/Predictor Set Intersection Analysis Manuscript Files/"

calculateAndWriteHallmarkPathways <- function(gmx,bg,genelist.jbsfl,jbsfl,genelist.dir){
  gs <- readLines(genelist.jbsfl)
  
  output <- lapply(gmx, function(x) {
    tab <- table(factor(bg %in% gs, levels = c(TRUE, FALSE)),
                 factor(bg %in% x, levels = c(TRUE, FALSE)))
    fit <- fisher.test(tab, alternative = "greater")
    interSection <- intersect(gs, x)
    interSection <- paste(interSection, collapse = ",")
    return(value = c(RATIO = as.numeric((diag(tab) / colSums(tab))[1]),
                     NOM_pval = fit$p.value,
                     INTERSECT = interSection))
  })
  
  fDF <- do.call(what = rbind, args = output) %>%
    as.data.frame() %>%
    mutate(NAME = names(gmx))
  
  fDF <- cbind(fDF,ADJ_pval=p.adjust(as.numeric(as.character(fDF$NOM_pval)), 
                            method = "BH"))
  # save results
  outFile <- fDF %>% 
    mutate(NEG_LOG10_FDR = -log10(ADJ_pval)) %>%
    arrange(NEG_LOG10_FDR) %>%
    dplyr::select(NAME, RATIO, NOM_pval, ADJ_pval, NEG_LOG10_FDR, INTERSECT)
  outFileName <- paste(genelist.dir,jbsfl,"HallmarkDatabasePathways.txt", sep="")
  write.table(outFile, 
              file = outFileName,
              quote = FALSE,
              sep = "\t",
              row.names = FALSE)
}

jbsfl1 <- "JBSFL1"
jbsfl2 <- "JBSFL2"
jbsfl3 <- "JBSFL3"
jbsfl4 <- "JBSFL4"
jbsfl5 <- "JBSFL5"
jbsfl6 <- "JBSFL6"
jbsfl7 <- "JBSFL7"
jbsfl8 <- "JBSFL8"
jbsfl9 <- "JBSFL9"
jbsfl10 <- "JBSFL10"
jbsfl11 <- "JBSFL11"
jbsfl12 <- "JBSFL12"

genelist.jbsfl1 <- paste(genelist.dir,jbsfl1,".txt",sep="")
genelist.jbsfl2 <- paste(genelist.dir,jbsfl2,".txt",sep="")
genelist.jbsfl3 <- paste(genelist.dir,jbsfl3,".txt",sep="")
genelist.jbsfl4 <- paste(genelist.dir,jbsfl4,".txt",sep="")
genelist.jbsfl5 <- paste(genelist.dir,jbsfl5,".txt",sep="")
genelist.jbsfl6 <- paste(genelist.dir,jbsfl6,".txt",sep="")
genelist.jbsfl7 <- paste(genelist.dir,jbsfl7,".txt",sep="")
genelist.jbsfl8 <- paste(genelist.dir,jbsfl8,".txt",sep="")
genelist.jbsfl9 <- paste(genelist.dir,jbsfl9,".txt",sep="")
genelist.jbsfl10 <- paste(genelist.dir,jbsfl10,".txt",sep="")
genelist.jbsfl11 <- paste(genelist.dir,jbsfl11,".txt",sep="")
genelist.jbsfl12 <- paste(genelist.dir,jbsfl12,".txt",sep="")

gmx <- "h.all.v6.0.symbols.gmt"
gmxFile <- paste(data.dir, gmx,sep="")
colNames <- max(count.fields(file = gmxFile, sep = "\t"))
colNames <- seq(from = 1, to = colNames)
colNames <- as.character(colNames)
gmx <- read.table(file      = gmxFile,
                  sep       = "\t",
                  quote     = "\"",
                  fill      = TRUE,
                  col.names = colNames,
                  row.names = 1)
gmx <- gmx[, -1]
gmx <- apply(gmx, MARGIN = 1, FUN = function(x) {
  return(value = setdiff(unname(x), ""))
})
names(gmx) <- toupper(names(gmx))

# extract the name of db from the gmx arguments
db <- basename(gmxFile)
db <- gsub(pattern = "[.].+$",
           replacement = "",
           db)
# obtain background
bg <- unique(unlist(gmx))
fisherIndex <- NULL

calculateAndWriteHallmarkPathways(gmx,bg,genelist.jbsfl1,jbsfl1,genelist.dir)
calculateAndWriteHallmarkPathways(gmx,bg,genelist.jbsfl2,jbsfl2,genelist.dir)
calculateAndWriteHallmarkPathways(gmx,bg,genelist.jbsfl3,jbsfl3,genelist.dir)
calculateAndWriteHallmarkPathways(gmx,bg,genelist.jbsfl4,jbsfl4,genelist.dir)
calculateAndWriteHallmarkPathways(gmx,bg,genelist.jbsfl5,jbsfl5,genelist.dir)
calculateAndWriteHallmarkPathways(gmx,bg,genelist.jbsfl6,jbsfl6,genelist.dir)
calculateAndWriteHallmarkPathways(gmx,bg,genelist.jbsfl7,jbsfl7,genelist.dir)
calculateAndWriteHallmarkPathways(gmx,bg,genelist.jbsfl8,jbsfl8,genelist.dir)
calculateAndWriteHallmarkPathways(gmx,bg,genelist.jbsfl9,jbsfl9,genelist.dir)
calculateAndWriteHallmarkPathways(gmx,bg,genelist.jbsfl10,jbsfl10,genelist.dir)
calculateAndWriteHallmarkPathways(gmx,bg,genelist.jbsfl11,jbsfl11,genelist.dir)
calculateAndWriteHallmarkPathways(gmx,bg,genelist.jbsfl12,jbsfl12,genelist.dir)