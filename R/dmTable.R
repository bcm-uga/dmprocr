#'dmTable
#'
#'Generate differential methylation data table from beta values. It basically extract all tumorous samples and controls. Then it computes the difference between each tumorous and the mean of control.
#'
#'@param data A matrix of row methylation data with TCGA sample ids as column names and probes ids as row names.
#'@param tumoral_ref a vector of ids corresponding to tumoral samples.
#'@param control_ref a vector of ids corresponding to control samples.
#'
#'@example examples/example-dmRandomDataset.R
#'@example examples/example-dmTable.R
#'
#'
#'@export
dmTable <- function(data, tumoral_ref, control_ref) {
  meanControl <- rowMeans(data[, control_ref], na.rm = TRUE)
  data    <- data[, tumoral_ref]
  AllDM = data - meanControl  
  return(AllDM)
}