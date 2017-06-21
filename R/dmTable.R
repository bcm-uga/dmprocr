#'dmTable
#'
#'Generate differential methylation data table from beta values. It basically extract all tumorous samples and controls. Then it computes the difference between each tumorous and the mean of control.
#'
#'@param data A matrix of row methylation data with TCGA sample ids as column names and probes ids as row names.
#'@param platform a dataframe with metadata types as columns names and probes ids as row names.
#'@param exp_group a dataframe of metadata about the samples. Rows are ids and columns information type.
#'@param TumoralRef a vector of ids corresponding to tumoral samples.
#'@param ControlRef a vector of ids corresponding to control samples.
#'@param NARM will delete empty rows in data and platform if TRUE.
#'
#'@example examples/example-data.R
#'@example examples/example-dmTable.R
#'
#'
#'@export

dmTable <- function(data, platform, exp_group, TumoralRef, ControlRef,  NARM = TRUE){
  #delete row all na
  if(NARM){
    ind       <- apply(data, 1, function(x) all(is.na(x)))
    data      <- data[!ind, ]
    platform  <- platform[!ind, ]
  }
  
  #compute the mean of Control
  meanControl         <- rowMeans(data[, ControlRef], na.rm = TRUE)
  
  #keep only tum samples in data, exp_group...
  data      <- data[, which(colnames(data) %in% TumoralRef)]
  exp_group <- exp_group[which(rownames(exp_group) %in% TumoralRef), ]
  
  
  #Perfom differential methylation
  list  <- lapply(colnames(data), function(x){
    data[, x] - meanControl
  })
  
  AllDM               <- do.call(cbind, list)
  AllDM               <- as.data.frame(AllDM)
  colnames(AllDM)     <- colnames(data)
  rownames(AllDM)     <- rownames(data)
  
  #delete list 'cause it's heavy
  rm(list)
  
  #check dimensions constraints
  dim(AllDM)
  dim(exp_group)
  dim(platform)
  
  #And return a list of object
  met_studyAllMean    <- list(data = AllDM, exp_group = exp_group, platform = platform)
  
  return(met_studyAllMean)
  
}