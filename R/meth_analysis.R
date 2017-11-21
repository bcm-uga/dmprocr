#'interpolate_gene
#'
#'Build a interpolated signal of differential methylation value for a gene.
#'
#'
#' @param vec A numeric vector specifying differential methylation signal.
#' @param probes_pos is a vector of probes position on the chromosome. 
#' @param tss the transcription start site on the same chromosome.
#' @param strand is the strand from which the gene is red. 
#' @param win is the width of the window on the chromosome in bp where the function will fetch probes position and differential methylation value, default is 5000.
#' @param slide is the maximum width slide you'll alow when comparing two curve, default is 0.
#' @param interp.by is resolution at which the function interpolate the probes signal, default is 20.
#'@export
interpolate_gene = function(vec, probes_pos, tss, strand, win, slide, interp.by) {
  xf         <- seq(tss - win - slide,  tss + win + slide, by = interp.by)
  if (sum(!is.na(vec))==0) {
    return(rep(NA, length(xf)))
  }
  xp         <- probes_pos[order(probes_pos)]
  yp         <- vec[order(probes_pos)]
  idx = !is.na(yp)
  xp = xp[idx]
  yp = yp[idx]
  # xp_orig = xp
  # yp_orig = yp
  #concatenate missing bp to x axis (in order to have 10000 bp windows everytime)
  xpA = xpA = c()
  if (tss - win - slide < min(xp)) {
    xpA = tss - win - slide
  }
  if (tss + win + slide > max(xp)) {
    xpB = tss + win + slide
  }
  # xpA <- seq(from = tss - win - slide, to = min(xp)-1 )
  # xpB <- seq(from = max(xp)+1, to = tss + win + slide)
  xp  <- c(xpA, xp, xpB)
  #add fictiv value to added bp in y axis
  yp  <- c( rep(0, length(xpA)), yp, rep(0, length(xpB)) )
  ######
  yf         <- signal::pchip(xp, yp, xf)
  # yf_orig = yf#         <- signal::pchip(xp_orig, yp_orig, xf)
  # layout(matrix(1:2, 1), respect=TRUE)
  # plot(xf, yf, ylim=range(c(yf_orig, yf)))
  # points(xp_orig, yp_orig, pch=16, col=2)
  # plot(xf, yf_orig, ylim=range(c(yf_orig, yf)))
  # points(xp_orig, yp_orig, pch=16, col=2)
  
  
  if (strand == "-") {
    yf = rev(yf)
  }
  return(yf)
}



#'compute_gene_meth_profile 
#'
#'Generate  differential methylation profile for a gene around the transcription start side 
#'
#' @param meth_data A matrix of row methylation data with TCGA sample ids as column names and probes ids as row names.
#' @param meth_platform a data.frame with metadata types as columns names and probes ids as row names.
#' @param gene A list that describe gene in bed format.
#' @param win is the width of the window on the chromosome in bp where the function will fetch probes position and differential methylation value
#' @param slide is the maximum width slide you'll alow when comparing two curve
#' @param interp.by is resolution at which the function interpolate the probes signal
#' @param pf_chr_colname string matching the name of the column in the platform that contain the chromosome information of probes
#' @param pf_pos_colname string matching the name of the column in the platform that contain the position information of probes
#' @param apply_func Function that will be used for apply.
#'@export
compute_gene_meth_profile = function(gene, meth_data, meth_platform, pf_pos_colname, pf_chr_colname, win, slide, interp.by, apply_func=apply) {  
  probe_idx = get_probe_names(gene   , 
    meth_platform=meth_platform      , 
    pf_pos_colname=pf_pos_colname    ,   
    pf_chr_colname=pf_chr_colname    , 
    up_str=win+slide                 , 
    dwn_str=win+slide 
  )
  
  if (length(probe_idx) == 0) {
    warning(paste0("No probes for gene ", gene[[4]],"(",gene[[5]],")."))
    return(NULL)
  } else {

    data         = meth_data[probe_idx,]
    probes_pos   = meth_platform[probe_idx, pf_pos_colname]
    strand       = gene[[6]]
    tss          = ifelse (strand == "+", as.numeric(gene[[2]]), as.numeric(gene[[3]]))

    # profile = dmProfile(gene_study_info, slide=slide, pf_pos_colname=pf_pos_colname)
    # return(profiles)
    if  (length(probe_idx) == 1) {
      data = t(data)
    }

    big2 = apply_func(
      data , 2,
      # vec = data[,5]
      interpolate_gene      ,
      probes_pos=probes_pos ,
      tss=tss               ,
      strand=strand         ,
      win=win               ,
      slide=slide           , 
      interp.by=interp.by
    )      


    meanOfInter2 = apply(big2, 1, mean, na.rm=TRUE)
    # plot(meanOfInter2, meanOfInter[,2])
    # meanOfInter2 == meanOfInter[,2]
    varOfInter2 = apply(big2, 1, var, na.rm=TRUE)
    # plot(varOfInter2, varOfInter[,2])
    # varOfInter2 == varOfInter[,2]

    xf         <- seq(tss - win - slide, tss + win + slide, by = interp.by)

    pond = as.numeric(sapply(xf, function(x) {
      # 0
      min(abs(x - probes_pos)) <= 20*interp.by
    }))
    pond[varOfInter2==0] = 0

    profile = data.frame(x=xf, y=meanOfInter2, var=varOfInter2, pond=pond, id=gene[[4]])
  
    return(profile)
  }
}

#' A Function That Extracts Probe Names of a Corresponding Gene from Platform Data
#'
#' @param gene A vector describing the gene (line of a bed file).
#' @param meth_platform A data frame describing CpG positions.
#' @param up_str   An integer specifying up stream size (in bp).
#' @param dwn_str  An integer specifying down stream size (in bp).
#' @param pf_pos_colname string matching the name of the column in the platform that contain the chromosome on which we find a probes.
#' @param pf_pos_colname string matching the name of the column in the platform that contain the position information of probes.
#' @return A vector of probe names
#' @export
get_probe_names = function(
  gene                                ,
  meth_platform                       , 
  pf_chr_colname="Chromosome"         ,
  pf_pos_colname="Start"              ,
  up_str=5000                         , 
  dwn_str=5000                        
 ) {
  # get gene properties
  chr =            gene[[1]]
  strand =         gene[[6]]
  gene_name =      gene[[4]]
  beg = as.numeric(gene[[2]])
  end = as.numeric(gene[[3]])
  
  # get meth infos
  if (strand == "-") {
    off_set_beg = dwn_str
    off_set_end = up_str
    tss = end
  } else {
    off_set_beg = up_str
    off_set_end = dwn_str
    tss = beg
  }

  ## Compute probes associated with the gene 
    probe_idx = rownames(meth_platform)[
      !is.na(meth_platform[[pf_pos_colname]]) & !is.na(meth_platform[[pf_chr_colname]]) &
      meth_platform[[pf_chr_colname]] == chr &
      meth_platform[[pf_pos_colname]] >= tss-up_str &
      meth_platform[[pf_pos_colname]] < tss+dwn_str
    ]    

  if (length(probe_idx) == 0) {
    warning(paste0("No probes for gene ", gene[[4]],"(",gene[[5]],")."))
    return(NULL)
  } else {
    return(probe_idx)    
  }
}

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
#'@example examples/example-dmDistance.R
#'
#'
#'@export
dmDistance <- function(dmprofileList, frechet = FALSE){
  m  = sapply(profiles, "[[","y")
  v  = sapply(profiles, "[[","var")
  k  = sapply(profiles, "[[","pond")
  d = matrix(0, nrow=ncol(m), ncol=ncol(m) )
  for (i in 1:(ncol(m)-1)) {
    for (j in (i+1):ncol(m)) {
      idx = k[,i] & k[,j]
      euc = (m[,i] - m[,j])^2 / (v[,i] + v[,j])
      # res = sum(euc * k[,i] * k[,j], na.rm=TRUE) / sum(k[,i] * k[,j], na.rm=TRUE)
      # res = sum(euc[idx] * k[idx,i] * k[idx,j]) / sum(k[idx,i] * k[idx,j])
      res = sum(euc[idx] * k[idx,i] * k[idx,j]) / sum(idx)
      res = sqrt(res)
      d[i,j] = res
      d[j,i] = res
    }
  }
  return(d)
}

