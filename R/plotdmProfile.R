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
    
    p1 <- ggplot2::ggplot(profile) +
      ggplot2::geom_line(ggplot2::aes_string(x = "x", y= "y"), size = 0.5, color = "red") +
      ggplot2::geom_ribbon(ggplot2::aes_string(x= "x", ymin= "y" - stats::sd, ymax= "y" + stats::sd), fill = "grey70", alpha = 0.8) +
      #theme
      ggplot2::geom_vline(xintercept = 0, alpha = 0.8) +
      ggplot2::geom_hline(yintercept = 0, alpha= 0.8) + 
      ggplot2::coord_cartesian(ylim = c(-1,1), xlim = c(min(profile$x), max(profile$x))) +
      ggplot2::theme(legend.position="none") +
      ggplot2::ggtitle(name)
    
    return(p1)
  })
  
  
}