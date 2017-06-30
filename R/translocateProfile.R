#'translocateProfile
#'
#'Produce a list of length bslide + 1 where each element is a dmProfile dataframe translocated around the tss 
#'
#'@param m1 is a dmProfile dataframe to translocate using bslide
#'@param bslide is the slide value divided by the width of each bin in dmProfile (= by.interp)
#'@param bwin is the window parameter divided by the width of each bin in dmProfile (= by.interp)
#'
#'@examples
#'#[warning]this function is not to be used outside of the context of dmDistance_translocate
#'
#'@export
translocateProfile <- function(m1, bslide, bwin){
  
  listmi <- list()
  
  for (k in -bslide:bslide){ 
    
    listmi[[k+bslide+1]]  <- m1[(1+bslide+k):(bslide+2*bwin+k), ]
    
  }
  
  return(listmi)
}