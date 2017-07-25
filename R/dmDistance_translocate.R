#'dmDistance_translocate
#'
#'Produce a list of two matrix : The distance matrix from a list of dmProfile. Each profile is translocated N times against another, we then keep the min(distance) in the matrix. The second matrix indicates which translocation returned the min(distance)
#'
#'@param dmprofileList a list of dmProfile
#'@param win is the width of the window on the chromosome in bp where the function will fetch probes position and differential methylation value, default is 5000.
#'@param slide is the maximum width slide you'll alow when comparing two curve, default is 0.
#'@param by.interp is resolution at which the function interpolate the probes signal, default is 20.
#'
#'@example examples/example-dmRandomDataset.R
#'@example examples/example-dmTable.R
#'@example examples/example-getalldmProfile.R
#'@example examples/example-dmDistance_translocate.R
#'
#'@export
dmDistance_translocate <- function(dmprofileList, win=5000, slide=500, by.interp = 20){
  
  #transform bp in bins of bp.length = by.interp
  bwin   <- win / by.interp
  bslide <- slide / by.interp 
  
  #create empty matrix
  m  <- matrix(rep(NA, length(dmprofileList)*length(dmprofileList)), nrow = length(dmprofileList), ncol = length(dmprofileList))
  mK <- matrix(rep(NA, length(dmprofileList)*length(dmprofileList)), nrow = length(dmprofileList), ncol = length(dmprofileList))
  
  for(i in 1:(length(dmprofileList))){ 
    
    m1 <- dmprofileList[[i]]
    
    listmi <- translocateProfile(m1, bslide, bwin) 
    
    # return(listmi)
    # for(j in (1+i):length(dmprofileList)){
    for(j in 1:length(dmprofileList)){
      
      m2 <- dmprofileList[[j]]
      m2 <- m2[(1+bslide):(bwin*2+bslide), ]
      
      #str(listmi)
      distvect <- sapply(listmi, eucli_dist, m2=m2)
      distvect <- sqrt(distvect)
      
      
      if(all(is.nan(distvect))){
        m[i,j]    <- NaN
        mK[i,j]   <- NaN
      }else{
        #distvect[k+slide+1]  <- sqrt(eucli.dist(m1, m2))
        m[i,j]    <- min(distvect, na.rm = TRUE)
        distvect[is.nan(distvect)] <- max(distvect, na.rm = TRUE)
        mK[i,j]   <- which.min(distvect) - bslide - 1
        # print(which.min(distvect))
      }
      
    }
  } 
  # dbmatrix <- list(dist=m, slide=mK)
  dbmatrix <-list(m = m, mK = mK)
  return(dbmatrix)
}
