#'getalldmProfile 
#'
#'This function is combination of getMethylInfo and dmProfile, it return a list of differential methylation profile. The bedline parameter is now a dataframe of bedline for which each will produce a dmProfile
#'
#'@param diff_meth_study list composed of a methylation differential table, an exp_grp dataframe and a platform 
#'@param bedfile a bedfile contaning information for one gene on each line
#'@param win is the width of the window on the chromosome in bp where the function will fetch probes position and differential methylation value
#'@param slide is the maximum width slide you'll alow when comparing two curve
#'@param interp.by is resolution at which the function interpolate the probes signal
#'@param pf_chr_colname string matching the name of the column in the platform that contain the chromosome on which we find a probes
#'@param pf_pos_colname string matching the name of the column in the platform that contain the position information of probes
#'
#'@example examples/example-dmRandomDataset.R
#'@example examples/example-dmTable.R
#'
#'@export
getalldmProfile <- function(diff_meth_study, bedfile, nbProbes = 0, win = 5000, slide = 0, interp.by = 20, pf_chr_colname="Chromosome", pf_pos_colname="Start") {
  
  genes <- bedfile$geneSymbol
  GenesInfo  <- lapply(1:nrow(bedfile), function(x) getMethylInfo(diff_meth_study = diff_meth_study, bedline = bedfile[x, ],pf_chr_colname =  pf_chr_colname, pf_pos_colname = pf_pos_colname))
  noprobes   <- which(sapply(GenesInfo, function(x){ is.data.frame(x$data) && nrow(x$data)==0 }))
  genes      <- genes[-noprobes]
  GenesInfo  <- GenesInfo[-noprobes] 
  
  
  toRemove <- sapply(GenesInfo, function(gene){
    
    if(nrow(gene$data) < nbProbes){
      TRUE
    } else{
      FALSE
    }
  })

  GenesInfo <- GenesInfo[!toRemove]
  genes <- genes[!toRemove]

  
  # allProfile <- lapply(GenesInfo, dmProfile, win=win, slide=slide, interp.by= interp.by)
  allProfile <- lapply(GenesInfo, function(resx){
    dmProfile(resx, slide=slide, pf_pos_colname = pf_pos_colname)
  }) 
  
}