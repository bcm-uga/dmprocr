#' getDifferentialTable
#'
#' Generate differential methylation data table from beta values. It basically extract all tumorous samples and controls. Then it computes the difference between each tumorous and the mean of control
#'
#' @param data A matrix of row methylation data with TCGA sample ids as column names and probes ids as row names
#' @param platform a dataframe with metadata types as columns names and probes ids as row names
#' @param exp_group a dataframe of metadata about the samples. Rows are ids and columns information type.
#'
#' @export  
#' @examples
#' library(gsubfn)
#' ###generate random samples ids  
#' regexTum <- "^TCGA-.{2}-.{4}-01.{1}-.{3}-.{4}-.{2}$"
#' regexWT <- "^TCGA-.{2}-.{4}-11.{1}-.{3}-.{4}-.{2}$"
#'
#' sampleId <- c(
#'          sapply(1:6, function(x) {
#'              gsub("\\^|\\$", "",
#'              gsubfn("\\.{([[:digit:]]+)}",
#'              ~ paste(rep(sample(c(LETTERS, 1:9), 1), n), collapse=""), regexTum)) }),
#'          sapply(1:4, function(x) {
#'              gsub("\\^|\\$", "",
#'              gsubfn("\\.{([[:digit:]]+)}", 
#'              ~ paste(rep(sample(c(LETTERS, 1:9), 1), n), collapse=""), regexWT)) })
#'              )
#' idProbes <- paste0(rep("probe", 10), 1:10) 
#' data <- matrix(runif(100,0,1), 10, 10)
#' colnames(data) <- sampleId
#' rownames(data) <- idProbes
#'
#' ###generate exp_group
#' exp_group <- as.data.frame(matrix(rep(1, 50),  10, 5), row.names = sampleId) 
#'
#' ###generate platfrom
#' platform  <-as.data.frame(matrix(rep(1, 50),  10, 5), row.names = idProbes)
#' 
#' ###get a methylation differential table
#' diff_table <- getDifferentialTable(data, platform, exp_group)



getDifferentialTable <- function(data, platform, exp_group){
  
  ## get the names of ind by type
  TumoralRef       <- grep("^TCGA-.{2}-.{4}-01.{1}-.{3}-.{4}-.{2}$", colnames(data), perl = TRUE  ,value = TRUE)
  ControlRef       <- grep("^TCGA-.{2}-.{4}-11.{1}-.{3}-.{4}-.{2}$", colnames(data), value = TRUE)
  OthersRef        <- grep("^TCGA-.{2}-.{4}-02", colnames(data), value = TRUE)
  
  #delete row all na
  ind                 <- apply(data, 1, function(x) all(is.na(x)))
  data      <- data[!ind, ]
  platform  <- platform[!ind, ]  
  
  
  #compute the mean of Control
  meanControl         <- rowMeans(data[, ControlRef], na.rm = TRUE)
  
  
  #keep only tum samples in data, exp_group...
  data      <- data[, which(colnames(data) %in% TumoralRef)]
  exp_group <- exp_group[which(rownames(exp_group) %in% TumoralRef), ] 
  platform  <- platform[!ind, ]
  
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
