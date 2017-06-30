#'eucli.dist 
#'
#'A function to compute adapted eucliddean distance on differential methylation profile. This function is called by dmDistance and dmDistance_translocate functions. 
#'
#'@param m1 a dmProfile to be compared with another one.
#'@param m2 the other one dmProfile.
#'
#'@example examples/example-dmRandomDataset.R
#'@example examples/example-dmTable.R
#'@example examples/example-getalldmProfile.R
#'@example examples/example-eucli_dist.R
#'
eucli.dist <- function(m1, m2){ #where m1 and m2 are methylation profile
  ##here we go, first the distance
  D.probe <-  (m1$y - m2$y)^2 / (m1$var + m2$var)
  
  ##then the moderation
  D.gene  <- sum(D.probe * m1$pond * m2$pond, na.rm = TRUE) / sum(m1$pond * m2$pond, na.rm = TRUE)
  
  return(D.gene)
  
}
