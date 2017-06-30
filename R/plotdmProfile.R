#'plotdmProfile
#'
#'Produce a list of differential methylation profile plot with ggplot2.
#'
#'@param dmprofileList a list of dmProfile 
#'
#'@example examples/example-dmRandomDataset.R
#'@example examples/example-dmTable.R
#'@example examples/example-getalldmProfile.R
#'@example examples/example-plotdmProfile
#'
#'
#'
#'@export
plotdmProfile <- function(dmprofileList){
  
  plots  <- lapply(dmprofileList, function(profile){  
    profile$sd <- sqrt(profile[, 3])
    name <- profile[, 5]
    
    p1 <- ggplot(profile) +
      geom_line(aes(x = x, y= y), size = 0.5, color = "red") +
      geom_ribbon(aes(x= x, ymin= y - sd, ymax= y + sd), fill = "grey70", alpha = 0.8) +
      #theme
      geom_vline(xintercept = 0, alpha = 0.8) +
      geom_hline(yintercept = 0, alpha= 0.8) + 
      coord_cartesian(ylim = c(-1,1), xlim = c(min(profile$x), max(profile$x))) +
      theme(legend.position="none") +
      ggtitle(name)
    
    return(p1)
  })
  
  
}