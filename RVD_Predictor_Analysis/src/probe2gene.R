################################################################################
## This code is suppose to convert a matrix of microarray probe measurements  ##
## to gene level data.                                                        ##
################################################################################
##                                                                            ##
## Author: Mehrad Mahmoudian                                                  ##
## Email: m.mahmoudian@gmail.com                                              ##
## License: MIT                                                               ##
##                                                                            ##
################################################################################
##                                                                            ##
## Copyright (c) 2017 Mehrad Mahmoudian                                       ##
##                                                                            ##
##     Permission is hereby granted, free of charge, to any person obtaining  ##
## a copy of this software and associated documentation files (the "Software" ##
## ), to deal in the Software without restriction, including without          ##
## limitation the rights to use, copy, modify, merge, publish, distribute,    ##
## sublicense, and/or sell copies of the Software, and to permit persons to   ##
## whom the Software is furnished to do so, subject to the following          ##
## conditions:                                                                ##
##                                                                            ##
##   The above copyright notice and this permission notice shall be included  ##
## in all copies or substantial portions of the Software.                     ##
##                                                                            ##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR ##
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   ##
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    ##
## THE, AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR      ##
## OTHER, LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,     ##
## ARISING FROM,, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     ##
## OTHER DEALINGS IN THE, SOFTWARE.                                           ##
##                                                                            ##
################################################################################

megre_probe2gene <- function(GE_df, GPL, method = "median", keep.one2many = F,
                             genesymbol_delimiter = " /// "){
    ## parameters:
    ##     GE_df:                the dataframe of gene expression (column names
    ##                               should be probe names)
    ##     GPL:                  the annotation file for probes. probe ID column
    ##                               should be names as "ID"
    ##     method:               the method which you want to choose columns
    ##                               (read probes) based on it
    ##     keep.one2many:        logical. indicate whether you want to keep
    ##                               those probes that map to more than one gene
    ##     genesymbol_delimiter: A character vector of length 1 to represent the
    ##                               delimiter that have been used in GPL to
    ##                               separate gene symbols in each row
    
    #-------[ loading the required packages ]-------#
    {
        require("progress")
        require("varhandle")
    }
    
    
    #-------[ pre-process ]-------#
    {
        # check if the symbol column exists in the given data
        if (!any(grepl(pattern = "symbol", x = colnames(GPL), ignore.case = T))) {
            stop(":v: There is no \"Symbol\" column")
        }else{
            colnames(GPL)[grepl(pattern = "symbol", x = colnames(GPL), ignore.case = T)] <- "Symbol"
        }
        
        
        ## make sure names are in correct compatible format
        GPL$ID <- make.names(GPL$ID)
        colnames(GE_df) <- make.names(colnames(GE_df))
        
        
        ## convert to dataframe
        if(class(GE_df) != "data.frame"){
            message(":!: Convert the GE_df to data.frame ...")
            GE_df <- data.frame(GE_df)
        }
        
        
        # remove those probes that does not exist in the GE but do exist in GPL
        GPL <- GPL[match(colnames(GE_df), GPL$ID), ]
        
        ## check if everything went right
        if (ncol(GE_df) != nrow(GPL)) {
            stop(":v: something went wrong while selecting the probes from GPL!")
        }
        if (any(colnames(GE_df) != GPL$ID)) {
            stop(":v: The order is not correct")
        }
        
        ## only keep the symbols we want
        # remove duplicates and NA
        GPL_symbol <- na.omit(unique(GPL$Symbol))
        # remove those that are empty ""
        GPL_symbol <- GPL_symbol[GPL_symbol != "" & GPL_symbol != "---"]
        
        # if user wants to keep the probes that map to many
        if(keep.one2many){
            # break the gene symbol list with more than one gene to multiple items
            GPL_symbol <- unique(unlist(strsplit(x = GPL_symbol, split = genesymbol_delimiter), use.names = F))
        }else{
            # remove those items that contains genesymbol_delimiter in their gene symbol list (meaning they are mapping to more than one gene)
            GPL_symbol <- GPL_symbol[-grep(x = GPL_symbol, pattern = genesymbol_delimiter)]
        }
        
        
        ## make sure there is no factor
        if(any(sapply(GE_df, class) == "factor")){
            message(":!: Unfactoring ...")
            GE_df <- unfactor(GE_df)
        }
        
    }
    
    
    #-------[ conversion ]-------#
    {
        # collective dataframe (to be filled in for loop)
        collective_df <- data.frame(row.names = row.names(GE_df))
        
        # inform user
        message(" :!: Processing the conversion ...")
        
        # initialize progressbar
        pb <- progress_bar$new(total = length(GPL_symbol))
        
        # iterate through the GPL_symbol
        for (g in GPL_symbol) {
            # create a subset for this specific gene
            tmp_probe_names <- GPL[grep(pattern = g, x = GPL$Symbol), "ID"]
            
            # if there is no probe that match to this gene!
            if (length(tmp_probe_names) == 0) {
                stop(":v: There is no probe that match to th gene", g, ", but there should be. something is wrong!")
            }
            
            # check the first character of each probe and it if is a number, add X infront of probe name
            tmp_probe_names[unlist(lapply(strsplit(x = tmp_probe_names, split = ""),
                                          function(x){check.numeric(x[1])}))] <- paste0("X", tmp_probe_names[unlist(lapply(strsplit(x = tmp_probe_names, split = ""),
                                                                                                                           function(x){check.numeric(x[1])}))])
            
            if (any(!is.element(tmp_probe_names, colnames(GE_df)))) {
                stop(":v: the probes we are looking for does not exist in the GE_df")
            }
            
            # create a subset from the GE_df
            tmp_subset <- GE_df[, tmp_probe_names]
            
            # if there were more than one probe for this particular gene
            if (length(tmp_probe_names) > 1) {
                # get the median or any method user chooses for each column
                tmp_median <- sapply(tmp_subset, function(x){
                    eval(call(method, x))
                })
                
                # overwrite the tmp_subset with the chosen column
                tmp_subset <- tmp_subset[, which.max(tmp_median)]
            }
            
            # append the new column to the collective dataframe
            collective_df <- cbind.data.frame(collective_df, tmp_subset)
            
            # update progressbar
            pb$tick()
        }
        
        
        # set column names of the collective dataframe to follow the gene symbols we have
        colnames(collective_df) <- GPL_symbol
        
        # inform user
        message(":^: Conversion was successful. The function got", ncol(GE_df), "probes as input and converted them to", ncol(collective_df), "genes.")
    }
    
    return(collective_df)
}


# #-------[ example ]-------#
# {
#     ## load packages
#     library("readr")
#     
#     # read data in
#     hgu133a_2 <- data.frame(read_csv("~/Downloads/HG-U133A_2.na35.annot.csv/HG-U133A_2.na35.annot.csv", comment = "#"))
#     expression_data <- data.frame(read_tsv("~/Downloads/ViralChallenge_training_EXPRESSION_RMA.tsv", comment = "#"))
#     
#     ## fix the dimention and probe names
#     row.names(expression_data) <- make.names(expression_data$FEATUREID)
#     expression_data <- expression_data[, -which(colnames(expression_data) == "FEATUREID")]
#     expression_data <- t.data.frame(expression_data)
#     
#     # fix the column name
#     colnames(hgu133a_2)[1] <- "ID"
#     
#     expression_data_gene_level <- megre_probe2gene(GE_df = expression_data, GPL = hgu133a_2)
# }


