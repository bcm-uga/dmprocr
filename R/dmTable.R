#'dmTable
#'
#'Generate differential methylation data table from beta values. It basically extract all tumorous samples and controls. Then it computes the difference between each tumorous and the mean of control.
#'
#'@param data A matrix of row methylation data with TCGA sample ids as column names and probes ids as row names.
#'@param platform a dataframe with metadata types as columns names and probes ids as row names.
#'@param exp_grp a dataframe of metadata about the samples. Rows are ids and columns information type.
#'@param tumoral_ref a vector of ids corresponding to tumoral samples.
#'@param control_ref a vector of ids corresponding to control samples.
#'@param RM_NA_ROWS will delete empty rows in data and platform if TRUE.
#'
#'@example examples/example-dmRandomDataset.R
#'@example examples/example-dmTable.R
#'
#'
#'@export
dmTable <- function(data, platform, exp_grp, tumoral_ref, control_ref,  RM_NA_ROWS=TRUE){
  #delete row all na
  if(RM_NA_ROWS){
    ind       <- apply(data, 1, function(x) all(is.na(x)))
    data      <- data[!ind, ]
    platform  <- platform[!ind, ]
  }
  
  #compute the mean of Control
  meanControl         <- rowMeans(data[, control_ref], na.rm = TRUE)
  
  #keep only tum samples in data, exp_grp...
  data    <- data[, colnames(data) %in% tumoral_ref]
  exp_grp <- exp_grp[colnames(data), ]
  
  
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
  dim(exp_grp)
  dim(platform)
  
  #And return a list of object
  diff_meth_study    <- list(data=AllDM, exp_grp=exp_grp, platform=platform)
  
  return(diff_meth_study)
  
}