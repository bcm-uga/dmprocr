#'getDifferentialTable
#'
#'This function generate a differential methylation data table from raw data. 
#'
#'@param data A matrix of row methylation data with TCGA sample ids as column names and probes ids as row names
#'@param platform a dataframe with metadata types as columns names and probes ids as row names
#'@param exp_group a dataframe of metadata about the samples. Rows are ids and columns information type.
#'
#'@export  


getDifferentialTable <- function(data, platform, exp_group){
  
  ## get the names of ind by type
  getTumoralRef       <- grep("^TCGA-.{2}-.{4}-01.{1}-.{3}-.{4}-.{2}$", colnames(data), perl = TRUE  ,value = TRUE)
  getControlRef       <- grep("^TCGA-.{2}-.{4}-11.{1}-.{3}-.{4}-.{2}$", colnames(data), value = TRUE)
  getOthersRef        <- grep("^TCGA-.{2}-.{4}-02", colnames(data), value = TRUE)
  
  #delete row all na
  ind                 <- apply(data, 1, function(x) all(is.na(x)))
  data      <- data[!ind, ]
  platform  <- platform[!ind, ]  
  
  
  #compute the mean of Control
  meanControl         <- rowMeans(data[, getControlRef], na.rm = TRUE)
  
  
  #delete control & unwanted tum samples from data, exp_group...
  data      <- data[, -which(colnames(data) %in% getControlRef)]
  exp_group <- exp_group[-which(rownames(exp_group) %in% getControlRef), ] 

  
  list  <- lapply(colnames(data), function(x){
    data[, x] - meanControl
  })
  
  AllDM               <- do.call(cbind, list)
  AllDM               <- as.data.frame(AllDM)
  colnames(AllDM)     <- colnames(data)
  rownames(AllDM)     <- rownames(data)
  
  #delete list 'cause it's heavy
  rm(list)
  
  
  dim(AllDM)
  dim(exp_group)
  dim(platform)
  
  
  met_studyAllMean    <- list(data = AllDM, exp_group = exp_group, platform = platform)
  
  return(met_studyAllMean)
  
}