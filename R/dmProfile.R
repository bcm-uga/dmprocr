#'dmProfile 
#'
#'Generate  differential methylation profile for genes around the transcription start side 
#'
#'@param gene_study_info is a subset of diff_meth_study list for a given gene bedline and a size of window  
#'@param win is the width of the window on the chromosome in bp where the function will fetch probes position and differential methylation value
#'@param slide is the maximum width slide you'll alow when comparing two curve
#'@param interp.by is resolution at which the function interpolate the probes signal
#'
#'@importFrom magrittr "%>%"
#'
#'
#'
#'@export
dmProfile <- function(gene_study_info, win = 5000, interp.by = 20, slide = 0){
  
  
  data        <- gene_study_info$data
  probes_start<- gene_study_info$probes$Start
  n           <- seq_along(gene_study_info$data)
  promoterPos <- gene_study_info$promoterPos
  strand      <- gene_study_info$strand
  
  
  #begin interpoling all 
  
  allsample <- lapply(n,
                      interpolateGene,
                      diff_table=data, probes_start=probes_start,
                      promoterPos=promoterPos,
                      strand=strand,
                      win=win,
                      slide=slide)
  big          <- do.call(rbind, allsample)
  
  
  meanOfInter <- stats::aggregate(big[, 2], list(big$position), mean)
  varOfInter  <- stats::aggregate(big[, 2], list(big$position), stats::var)
  colnames(meanOfInter) <- c("x", "y")
  meanOfInter$var <- varOfInter[, 2]
  
  
  #########build the density vector
  
  #get the position of probes again
  if(strand == "-"){
    probes_start <- 0 + (promoterPos - probes_start)
  } else {
    probes_start <- 0 - (promoterPos - probes_start)
  }
  
  # probes_start <-probes_start / 2500
  x <- meanOfInter$x
  
  #vector of interpolation + position of probes
  xxold <- sort(c(probes_start, x))
  
  #find index of old in bigvect
  indx <- match(probes_start, xxold) 
  
  #find index of point +-20 around a probes position 
  indxAround <- lapply(indx,  function(x) x + -20:20) %>% unlist() %>% intersect(seq_along(xxold))
  
  #initialized empty vector
  A <- double(length(xxold))
  
  #Give value to empty A vector
  A[indxAround] <- 1
  
  #delete xold from x2 with the indx
  A <- A[-indx]
  
  A[which(meanOfInter$y == "0")] <- 0 
  
  meanOfInter$pond <- A
  meanOfInter$id   <- rep(gene_study_info$id, length(meanOfInter$pond)) 
  
  return(meanOfInter)
  
}  