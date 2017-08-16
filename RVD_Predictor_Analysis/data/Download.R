# ---
# title: Download
# author: Joshua Burkhart
# date: Aug 10, 2017
# ---

# global variables
SYNAPSE_PARENT_ID <- "syn7067265"
SYNAPSE_QUERY <- paste("select * from file where file.parentId=='",SYNAPSE_PARENT_ID,"'",sep="")
DOWNLOAD_DIR <- "../data/downloads/"
MD5HASH_PATH <- paste(DOWNLOAD_DIR,"md5hash.txt",sep="")
CFG_PATH <- "../config/SynapseCredentials.dat"

if(!file.exists(MD5HASH_PATH)){
  
  # load library
  library(synapseClient)
  library(tools)
  
  # load credentials
  credentials <- read.csv(CFG_PATH,stringsAsFactors = FALSE)
  
  # login to Synapse
  synapseClient::synapseLogin(username = credentials$username, password =  credentials$password)
  
  # retrieve synapse IDs of all files of the parent synapse ID
  fileDF <- synapseClient::synQuery(SYNAPSE_QUERY)
  
  # download files
  lapply(fileDF$file.id, 
         FUN = synapseClient::synGet, 
         downloadLocation = DOWNLOAD_DIR)
  
  # get downloaded filenames
  downloaded_csvs <- paste(DOWNLOAD_DIR,list.files(DOWNLOAD_DIR,pattern = "*.csv"),sep="")
  
  # calculate md5sum
  md5hash <- tools::md5sum(downloaded_csvs)
  
  # store md5sum
  write(md5hash,file = MD5HASH_PATH)
  
}