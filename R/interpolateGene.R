#'InterpolateGene
#'
#'build a curve from differential methylation matrix 
#'
#'
#'@param n is a vector of index of the column of the data matrix (each corresponding to one sample).
#'@param diff_table is a differential methylation data matrix.
#'@param probes_start is a vector of probes position on the chromosome. 
#'@param promoterPos the transcription start site on the same chromosome.
#'@param strand is the strand from which the gene is red. 
#'@param win is the width of the window on the chromosome in bp where the function will fetch probes position and differential methylation value, default is 5000.
#'@param slide is the maximum width slide you'll alow when comparing two curve, default is 0.
#'@param interp.by is resolution at which the function interpolate the probes signal, default is 20.
#'
#'@export
interpolateGene <- function(n ,diff_table, probes_start, promoterPos, strand, win=5000, slide=0, interp.by=20){
  
  
  indx       <- diff_table[n]
  xp         <- probes_start[order(probes_start)]
  yp         <- indx[order(probes_start), ]
  
  
  #concatenate missing bp to x axis (in order to have 10000 bp windows everytime)
  xpA <- seq(from = promoterPos - win - slide, to = min(xp)-1 )
  xpB <- seq(from = max(xp)+1, to = promoterPos + win + slide)
  xp  <- c(xpA, xp, xpB)
  
  
  
  #add fictiv value to added bp in y axis
  yp  <- c( rep(0, length(xpA)), yp, rep(0, length(xpB)) )
  
  
  
  ######
  xf         <- seq(promoterPos - win - slide, promoterPos + win + slide, by = interp.by)
  yf         <- signal::pchip(xp, yp, xf)
  
  
  
  if(strand == "-"){
    xf <- 0 + (promoterPos - xf)
  } else {
    xf <- 0 - (promoterPos - xf)  
  }  
  
  # xf <- xf/2500
  
  bigStacslideed <- data.frame(position = xf, valuesInter = yf, ind = rep(names(indx), length(xf)))
  
  
}