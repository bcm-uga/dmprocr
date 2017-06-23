#'getMethylInfo
#'
#'Return a subset of diff_meth_study list for a given gene bedline and a size of window 
#'
#'@param diff_meth_study list composed of a methylation differential table, an exp_grp dataframe and a platform 
#'@param bedline is a line of of the bedfile dataframe
#'@param  win is the width of the window on the chromosome in bp where the function will fetch probes position and differential methylation value. Default is 5000 bp 
#'
#'@export
getMethylInfo <- function(diff_meth_study, bedline, win = 5000) {
  
  strand <- bedline[[5]]
  
  if(strand == "-"){
    txstart <- bedline[[4]]
  } else{
    txstart <- bedline[[3]]   
  }
  
  chrom      <- bedline[[2]]
  
  probes     <- diff_meth_study$platform[diff_meth_study$platform$Chromosome == chrom, ] 
  
  probes     <- diff_meth_study$platform[diff_meth_study$platform$Start < txstart + win & diff_meth_study$platform$Start > txstart - win, ]
  
  data       <-  diff_meth_study$data[which(rownames(diff_meth_study$data) %in% probes$Composite.Element.REF), ]
  data       <-  data[, !colSums(is.na(data)) > 0]
  gene_study_info   <- list(data = data, probes = probes, promoterPos = txstart, strand = strand, id = levels(droplevels(bedline[[1]])) )
  
  return(gene_study_info)
  
}