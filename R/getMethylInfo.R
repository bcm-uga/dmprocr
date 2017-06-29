#'getMethylInfo
#'
#'Return a subset of diff_meth_study list for a given gene bedline and a size of window
#'
#'@param diff_meth_study list composed of a methylation differential table, an exp_grp dataframe and a platform
#'@param bedline is a line of of the bedfile dataframe
#'@param win is the width of the window on the chromosome in bp where the function will fetch probes position and differential methylation value. Default is 5000 bp
#'@param pf_chr_colname string matching the name of the column in the platform that contain the chromosome on which we find a probes
#'@param pf_pos_colname string matching the name of the column in the platform that contain the position information of probes
#'
#'@example examples/example-dmRandomDataset.R
#'@example examples/example-dmTable.R
#'@example examples/example-getMethylInfo.R
#'
#'@export
getMethylInfo <- function(diff_meth_study, bedline, win = 5000, pf_chr_colname="Chromosome", pf_pos_colname="Start") {
  
  strand <- bedline[[6]]
  if(strand == "-"){
    txstart <- as.numeric(bedline[[3]])
  } else{
    txstart <- as.numeric(bedline[[2]])   
  }
  
  chrom      <- bedline[[1]]
  idx = diff_meth_study$platform[[pf_chr_colname]] == chrom & diff_meth_study$platform[[pf_pos_colname]] < txstart + 5000 & diff_meth_study$platform[[pf_pos_colname]] > txstart - 5000
  probes = rownames(diff_meth_study$platform)[idx]
  data       <-  diff_meth_study$data[idx, ]
  data       <-  data[, !colSums(is.na(data)) > 0]
  gene_study_info   <- list(data = data, probes = diff_meth_study$platform[probes,], promoterPos = txstart, strand = strand, id = bedline[[4]])
  
  return(gene_study_info)
}