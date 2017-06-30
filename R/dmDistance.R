#'dmDistance
#'
#'Perform either a frechet distance from the kmlShape package or the eucliddean distance modified from this package (default) between two dmProfile. Return a distance matrix.
#'[warning]Carefully use the frechet distance as it can be heavy computing when dealing with large set of profile. Complexity of the profile also weight on the memory usage. 
#'
#'@param dmprofileList a list of dmProfile 
#'@param frechet a boolean specify if frechet distance will be computed.
#'
#'@example examples/example-dmRandomDataset.R
#'@example examples/example-dmTable.R
#'@example examples/example-getalldmProfile.R
#'@example examples/example-dmDistance
#'
#'
#'@export
dmDistance <- function(dmprofileList, frechet = FALSE){
  
  #create empty matrix
  m <- matrix(nrow = length(dmprofileList), ncol = length(dmprofileList))
  
  for(i in 1:length(dmprofileList)){ 
    
    for(j in 1:length(dmprofileList)){
      m1 <- dmprofileList[[i]]  
      m2 <- dmprofileList[[j]]
      
      
      if(frechet){
        
        frech.dist <- kmlShape::distFrechet(Px = m1$x, Py = m1$y, Qx = m2$x, Qy = m2$y)          
        m[i,j] <- frech.dist
        
      }else {
        
        m[i,j] <- sqrt(eucli_dist(m1, m2))
        
      }
      
    }     
  } 
  return(m)
}