# Authors: Magali Richard, CNRS
# magali.richard@univ-grenoble-alpes.fr
#
#---------------------------------------------
#' Perform differential analysis of RNAseq data.
#' 
#' Perform differential analysis of transcriptomic data (RNAseq) using DESeq2 R package.
#' @param data_trscr A \code{data} matrix that contains transcriptome information (RNAseq counts from HTseq). Columns correspond to indivuals, row correspond to genes.
#' @param exp_grp A \code{exp_grp} dataframe that contains metadatas on \code{data_trscr} individuals.
#' @param gene_list A \code{gene_list} bedfile containing the genes to screen for differential expression.
#' @param filter_indiv A vector of individual names to be screened for differential expression. Optionnal (set on "no_filter" by default).
#' @param alpha A parameter to indicate the significance cutoff used by \code{DESeq2::results} funtion for optimizing the independent filtering (by default 0.05). If the adjusted p-value cutoff (FDR) will be a value other than 0.1, alpha should be set to that value.
#' @param contrast A vector containing the constrast to be used to estimate the logarithmic fols change. By default: c("tissue_status","patho","normal")
#' @param fitType A DESeq fitting paramater, by default set to "parametric"
#' @param normalization_factor A matrix of normalization to preempt DESeq2 sizeFactors. Optionnal (set on "no_factor" by default).
#' 
#' @return A \code{gene_list} table including log2FoldChange and adjusted p-value (padj) computed by DESeq2 and a \code{data_ntrscr} matrix of normalized counts.
#' 
#' 
#' @importFrom DESeq2 "DESeqDataSetFromMatrix"
#' @importFrom DESeq2 "DESeq"
#' @importFrom DESeq2 "results"
#' @importFrom DESeq2 "counts"
#' @importFrom stats "as.formula"
#' @importFrom DESeq2 "normalizationFactors"
#'
#' @export

RNAseq_diffAnalysis = function(data_trscr, exp_grp, gene_list, filter_indiv = "no_filter", alpha = 0.05, contrast=c("tissue_status","patho","normal"), fitType="parametric", normalization_factor ="no_factor") {  
  if (filter_indiv[1] == "no_filter") {
    print("all individuals will be used in the differential analysis")
    filter_indiv = colnames(data_trscr)
  }
  if (length(which(filter_indiv %in% colnames(data_trscr))) != length(filter_indiv) ) {
    stop("ERROR, filter_indiv does not fit to indivuals present in the data matrix")
  }  
  #Generate the dds for DESeq2 (we start from a count matrix)
  countData <- data_trscr[rownames(gene_list),  filter_indiv] #data matrix
  colData = exp_grp[filter_indiv, ]
  #Check wether columns ids of the count matrix corresponds to rows ids of the column data 
  #if (length(which(rownames(colData) == colnames(countData))) != ncol(countData) ) {
  #  stop("ERROR, exp_grp does not correspond to data matrix")
  #  }
  #Construct a DESeqDataSet
  colData[[contrast[1]]] = as.factor(colData[[contrast[1]]])
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = countData,
                                colData = colData,
                                design = as.formula(paste("~", as.factor(contrast[1]))))
  
  # If required, preempt normalization step
  if ( normalization_factor == "no_factor") {
    print("No preemption of sizeFactors")
  } else {
    print("Preempt sizeFactors")
    DESeq2::normalizationFactors(dds) <- normalization_factor
  }
  
  #Perform the differential analysis 
  dds <- DESeq2::DESeq(dds, fitType= fitType)
  #Extract the results of differential analysis
  res <- DESeq2::results(dds,alpha = alpha ,  contrast=contrast)
  result_list = gene_list
  result_list$log2FoldChange = res[rownames(result_list), "log2FoldChange"]
  result_list$padj = res[rownames(result_list), "padj"]
  #Generate the normalized data matrix
  data_ntrscr = DESeq2::counts(dds, normalized=TRUE)
  return(list(result_list = result_list, data_ntrscr = data_ntrscr))
}

# Authors: Magali Richard, CNRS
# magali.richard@univ-grenoble-alpes.fr
#
#---------------------------------------------
#' Normalize RNAseq data according to DESeq2 method.
#' 
#' Perform normalization of transcriptomic data (RNAseq) using DESeq2 R package.
#' @param data_trscr A \code{data} matrix that contains transcriptome information (RNAseq counts from HTseq). Columns correspond to indivuals, row correspond to genes.
#' @param exp_grp A \code{exp_grp} dataframe that contains metadatas on \code{data_trscr} individuals.
#' @param gene_list A \code{gene_list} bedfile containing the genes to screen for differential expression.
#' @param filter_indiv A vector of individual names to be screened for differential expression. Optionnal (set on "no_filter" by default).
#' @param alpha A parameter to indicate the significance cutoff used by \code{DESeq2::results} funtion for optimizing the independent filtering (by default 0.05). If the adjusted p-value cutoff (FDR) will be a value other than 0.1, alpha should be set to that value.
#' @param contrast A vector containing the constrast to be used to estimate the logarithmic fols change. By default: c("tissue_status","patho","normal")
#' @param fitType A DESeq fitting paramater, by default set to "parametric"
#' 
#' @return A list composed of a \code{data_ntrscr} matrix of normalized counts and a \code{size_factor} vector of normalization factors.
#' 
#' 
#' @importFrom DESeq2 "DESeqDataSetFromMatrix"
#' @importFrom DESeq2 "DESeq"
#' @importFrom DESeq2 "results"
#' @importFrom DESeq2 "counts"
#' @importFrom stats "as.formula"
#'
#' @export

RNAseq_normalization = function(data_trscr, exp_grp, gene_list, filter_indiv = "no_filter", alpha = 0.05, contrast=c("tissue_status","patho","normal"), fitType="parametric") {  
  if (filter_indiv[1] == "no_filter") {
    print("all individuals will be used in the normalization process")
    filter_indiv = colnames(data_trscr)
  }
  if (length(which(filter_indiv %in% colnames(data_trscr))) != length(filter_indiv) ) {
    stop("ERROR, filter_indiv does not fit to indivuals present in the data matrix")
  }  
  #Generate the dds for DESeq2 (we start from a count matrix)
  countData <- data_trscr[rownames(gene_list),  filter_indiv] #data matrix
  colData = exp_grp[filter_indiv, ]
  #Check wether columns ids of the count matrix corresponds to rows ids of the column data 
  #if (length(which(rownames(colData) == colnames(countData))) != ncol(countData) ) {
  #  stop("ERROR, exp_grp does not correspond to data matrix")
  #  }
  #Construct a DESeqDataSet
  colData[[contrast[1]]] = as.factor(colData[[contrast[1]]])
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = countData,
                                        colData = colData,
                                        design = as.formula(paste("~", as.factor(contrast[1]))))

  #Generate the normalized data matrix
  dds = DESeq2::estimateSizeFactors(dds)
  data_ntrscr = DESeq2::counts(dds, normalized=TRUE)
  size_factor = DESeq2::sizeFactors(dds)
   
  return(list(data_ntrscr = data_ntrscr, size_factor = size_factor))
}


# Authors: Florent Chuffart, INSERM
#
#---------------------------------------------
#' Perform a Quick differential analysis of RNAseq data.
#' 
#' Perform a quick differential analysis of transcriptomic data (RNAseq) using non parametric test.
#' @param data_trscr A \code{data} matrix that contains transcriptome information (RNAseq counts from HTseq). Columns correspond to indivuals, row correspond to genes.
#' @param exp_grp A \code{exp_grp} dataframe that contains metadatas on \code{data_trscr} individuals.
#' @param gene_list A \code{gene_list} bedfile containing the genes to screen for differential expression.
#' @param filter_indiv A vector of individual names to be screened for differential expression. Optionnal (set on "no_filter" by default).
#' @param contrast A vector containing the constrast to be used to estimate the logarithmic fols change. By default: c("tissue_status","patho","normal")
#' @param apply_func Function to be used for apply.
#' 
#' @return A \code{gene_list} table including log2FoldChange and adjusted p-value (padj) computed by DESeq2 and a \code{data_ntrscr} matrix of normalized counts.
#' 
#' @importFrom DESeq2 "DESeqDataSetFromMatrix"
#' @importFrom DESeq2 "DESeq"
#' @importFrom DESeq2 "results"
#' @importFrom DESeq2 "counts"
#' @importFrom stats "as.formula"
#' @importFrom stats wilcox.test
#' @importFrom stats p.adjust
#'
#' @export

fad_RNAseq_diffAnalysis = function(data_trscr, exp_grp, gene_list, filter_indiv = "no_filter", contrast=c("tissue_status","patho","normal"), apply_func=apply) {  
  if (filter_indiv[1] == "no_filter") {
    print("all individuals will be used in the differential analysis")
    filter_indiv = colnames(data_trscr)
  }
  if (length(which(filter_indiv %in% colnames(data_trscr))) != length(filter_indiv) ) {
    stop("ERROR, filter_indiv does not fit to indivuals present in the data matrix")
  }  
  #Generate the dds for DESeq2 (we start from a count matrix)
  countData <- data_trscr[rownames(gene_list),  filter_indiv] #data matrix
  colData = exp_grp[filter_indiv, ]
  #Check wether columns ids of the count matrix corresponds to rows ids of the column data 
  #if (length(which(rownames(colData) == colnames(countData))) != ncol(countData) ) {
  #  stop("ERROR, exp_grp does not correspond to data matrix")
  #  }
  #Construct a DESeqDataSet
  colData[[contrast[1]]] = as.factor(colData[[contrast[1]]])
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = countData,
                                colData = colData,
                                design = as.formula(paste("~", as.factor(contrast[1]))))
  #Generate the normalized data matrix
  dds = DESeq2::estimateSizeFactors(dds)
  #Generate the normalized data matrix
  data_ntrscr = DESeq2::counts(dds, normalized=TRUE)

  idx1 = rownames(exp_grp)[!is.na(exp_grp[[contrast[1]]]) & exp_grp[[contrast[1]]] == contrast[2]]
  idx2 = rownames(exp_grp)[!is.na(exp_grp[[contrast[1]]]) & exp_grp[[contrast[1]]] == contrast[3]]
  idx1 = idx1[idx1 %in% colnames(data_ntrscr)]
  idx2 = idx2[idx2 %in% colnames(data_ntrscr)]
  w_pval = apply_func(data_ntrscr, 1, function(l) {
    # l = data_ntrscr[100,]
    wilcox.test(l[idx1], l[idx2])$p.value
  })
  # plot(density(log2(data_ntrscr)))
  m1 = apply(data_ntrscr[,idx1], 1, mean)
  m2 = apply(data_ntrscr[,idx2], 1, mean)  
  l2fc = log2(m1+1) - log2(m2+1)
  result_list = gene_list
  result_list$log2FoldChange = l2fc[rownames(result_list)] #res[rownames(result_list), "log2FoldChange"]
  result_list$pval = w_pval[rownames(result_list)]         #res[rownames(result_list), "padj"]
  result_list$padj = p.adjust(result_list$pval, method="BH")
  return(list(result_list = result_list, data_ntrscr = data_ntrscr))
}

#---------------------------------------------
#' Perform linear regression on expression and CNV data.
#' 
#' Perform linear regression on expression (counts) and CNV data.
#' @param data_ntrscr A \code{data} matrix that contains normalized RNAseq counts (from DESeq2 analysis). Columns correspond to indivuals, row correspond to genes.
#' @param data_cnv A \code{data} matrix that contains CNV data
#' @param exp_grp A \code{exp_grp} data.frame that contains metadatas on \code{data_trscr} individuals.
#' @param gene_list A \code{gene_list} bedfile containing the genes for which the linear regression will be perform.
#' @param filter_indiv A vector of individual names to be screened for differential expression. Optionnal, all individual if missing.
#' @param contrast A vector containing the constrast to be used to read metadata
#' @param cnv_filter A vector of two values indicating between which quantiles (for cnv data) the regression should be performed, by default set to c(0.025, 0.975). Should be set to FALSE if no filter is required.
#' @param apply_func Function to be used for apply.
#' 
#' @return A \code{gene_list} table including a z_score value associated with the linear model to the gene_list used in entry.
#'
#' 
#' @importFrom stats "quantile"
#' @importFrom stats "lm"
#'
#' @export



RNAseq_cnv_reg <- function(data_ntrscr, data_cnv, exp_grp, gene_list, filter_indiv, contrast=c("tissue_status","patho","normal"), cnv_filter = c(0.025, 0.975), apply_func=apply){
  if (missing(filter_indiv)) {
    print("all individuals will be used in the differential analysis")
    filter_indiv = exp_grp$dmprocr_ID
  }
  # Select tumor samples for individual that display cnv and trscr data
  cur = exp_grp[exp_grp[[contrast[1]]] == contrast[2], ] # selection of pathogical samples
  tc_indivs = rownames(cur)[cur$trscr == 1 & cur$cnv == 1 & cur$dmprocr_ID %in% filter_indiv]
  # Generate data.frame for linear regression with a z_score corresponding to the standardized beta coefficient 
  if (length(cnv_filter) == 1) {
    print("patient are not filtered upon cnv values")
  } else {
    print("patient are filtered upon cnv values")
  }
  
  data_cnv_transcr = cbind(data_ntrscr[rownames(gene_list), tc_indivs], data_cnv[rownames(gene_list), tc_indivs])
  z_scores = apply_func(data_cnv_transcr, 1, function (l) {
    cur_data = data.frame(patient=tc_indivs, trscr=l[1:length(tc_indivs)], cnv=l[(length(tc_indivs)+1):(2*length(tc_indivs))])

    # Check if values are different of ?? 0??
    if (sum(abs(cur_data$cnv)) == 0 ) {
      z_score = NA
    } else {
      #If no filter is required, perform the regression on all data
      if (length(cnv_filter) == 1) {
        red_cur_data = cur_data
      #If filter required, perform the regression on reduced data.frame
      } else {
        red_cur_data = cur_data[cur_data$cnv < quantile(cur_data$cnv, max(cnv_filter)) &
                                cur_data$cnv > quantile(cur_data$cnv, min(cnv_filter))  ,]
      }
      mod = lm(trscr ~ cnv, data=red_cur_data)
      #Extract the z_score from the model
      z_score = summary(mod)$coefficients[6]
    }

    # gene_id = gene
    return(z_score)
  })



  # z_scores = sapply(rownames(gene_list), function (gene) {
  #   # i = which(rownames(gene_list) == gene)
  #   # if (i %% 100 == 0){print (paste(gene,i, "/", length(rownames(gene_list)), sep=" "))}
  #   #Generate a temporary data.frame with counts values (trscr) and cnv value for each patient
  #   # tmp_data = apply (t(tc_indivs), 2, function (patient) {
  #   #   trscr = data_ntrscr[gene, patient]
  #   #   cnv = data_cnv[gene, patient]
  #   #   return(list(patient = patient, trscr=trscr, cnv = cnv))
  #   # })
  #   # cur_data = (do.call(rbind, tmp_data))
  #   # cur_data = data.frame(lapply(data.frame(cur_data, stringsAsFactors=FALSE), unlist), stringsAsFactors=FALSE)
  #   cur_data = data.frame(patient=tc_indivs, trscr=data_ntrscr[gene, tc_indivs], cnv=data_cnv[gene, tc_indivs])
  #
  #   # Check if values are different of ?? 0??
  #   if (sum(abs(cur_data$cnv)) == 0 ) {
  #     z_score = NA
  #   } else {
  #     #If no filter is required, perform the regression on all data
  #     if (length(cnv_filter) == 1) {
  #       red_cur_data = cur_data
  #     #If filter required, perform the regression on reduced data.frame
  #     } else {
  #       red_cur_data = cur_data[cur_data$cnv < quantile(cur_data$cnv, max(cnv_filter)) &
  #                               cur_data$cnv > quantile(cur_data$cnv, min(cnv_filter))  ,]
  #     }
  #     mod = lm(trscr ~ cnv, data=red_cur_data)
  #     #Extract the z_score from the model
  #     z_score = summary(mod)$coefficients[6]
  #   }
  #
  #   gene_id = gene
  #   return(z_score)
  # })
  
  #merging with existing gene_list
  if (length(z_scores) != nrow(gene_list) | sum(names(z_scores) != rownames(gene_list))!=0) {
    stop("ERROR, some genes are missing in the linear regression z_score")
  }
    
  result_list = gene_list
  result_list$z_score = z_scores[rownames(result_list)]
  return(result_list)
}

#---------------------------------------------
#' Calculate differential expression on a subset of patient.
#' 
#' Calculate differential expression on a subset of patient, in particular those with no CNV
#' @param data_ntrscr A \code{data} matrix that contains normalized RNAseq counts (from DESeq2 analysis). Columns correspond to indivuals, row correspond to genes.
#' @param data_cnv A \code{data} matrix that contains CNV data
#' @param exp_grp A \code{exp_grp} dataframe that contains metadatas on \code{data_trscr} individuals.
#' @param gene_list A \code{gene_list} bedfile containing the genes to screen for differential expression.
#' @param filter_indiv A vector of individual names to be screened for differential expression. Optionnal (set on "no_filter" by default).
#' @param contrast A vector containing the constrast to be used to estimate the logarithmic fols change. By default: c("tissue_status","patho","normal")
#' @param no_cnv_filter A vector of two values indicating the CNV values threshold for no CNV, by default set to c(-0.2, 0.2). 
#' 
#' 
#' @return A \code{gene_list} table including a noCNVlog2FC value associated with differential expression between sample without CNV.
#'
#' @export

noCNV_diffAnalysis = function(data_ntrscr, data_cnv, exp_grp, gene_list, filter_indiv = "no_filter", contrast=c("tissue_status","patho","normal"), no_cnv_filter = c(-0.2, 0.2)){
  noCNV_data = t(sapply(rownames(gene_list), function (gene) {
    # print(gene)
    # i = which(rownames(gene_list) == gene)
    # if (i %% 500 == 0){print (paste(gene,i, "/", length(rownames(gene_list)), sep=" "))}
    
    idx_nocnv = colnames(data_cnv)[data_cnv[gene,] > min(no_cnv_filter) 
                                 & data_cnv[gene,] < max(no_cnv_filter)]
    idx_trscr = rownames(exp_grp)[exp_grp$trscr == 1]
    cur = exp_grp[intersect(idx_trscr, idx_nocnv), ]
    idx_p = rownames(cur)[cur[[contrast[1]]] == contrast[2]] #look for pathological samples
    idx_c = rownames(cur)[cur[[contrast[1]]] == contrast[3]] #look for ctrl samples
    ntrscr_p = mean(as.numeric(data_ntrscr[gene,idx_p]))
    ntrscr_c = mean(as.numeric(data_ntrscr[gene,idx_c]))
    noCNVlog2FC = log2(ntrscr_p+1) - log2(ntrscr_c+1)
    return(c(noCNVlog2FC=noCNVlog2FC, nb_p=length(idx_p), nb_c=length(idx_c)))
  }))
  result_list = cbind(gene_list, noCNV_data[rownames(gene_list),])
  return(result_list)
}

#---------------------------------------------
#' Generate a candidate table
#' 
#' @param gene_list A \code{gene_list} bedfile containing the log2FC, zscore and noCNVlog2FC values.
#' @param padj_thresh A threshold p-value under which expression is significantly different by DESeq analysis
#' @param z_score_thresh A threshold zscore value, candidate genes should display a zcore value such as -t_zscore<zscore<t_zscore
#' @param noCNVlog2FC_thresh A threshold noCNVlog2FC, candidate genes should display a noCNVlog2FC value such as noCNVlog2FC > noCNVlog2FC_thresh or noCNVlog2FC < noCNVlog2FC_thresh
#' 
#' 
#' @return A \code{gene_list} table including of candidate genes
#' @export


select_candidates = function(gene_list, padj_thresh = 0, z_score_thresh, noCNVlog2FC_thresh){
  if (padj_thresh == 0) {
  print('selection with no pvalue')
  idx_z_score     = which(!is.na(gene_list$z_score    ) & abs(gene_list$z_score    ) < z_score_thresh)    
  idx_noCNVlog2FC = which(!is.na(gene_list$noCNVlog2FC) & abs(gene_list$noCNVlog2FC) > noCNVlog2FC_thresh)
  idx = intersect(idx_z_score, idx_noCNVlog2FC)
  candidate_list = gene_list[idx, ]
  } else {
  idx_z_score     = which(!is.na(gene_list$z_score    ) & abs(gene_list$z_score    ) < z_score_thresh)    
  idx_noCNVlog2FC = which(!is.na(gene_list$noCNVlog2FC) & abs(gene_list$noCNVlog2FC) > noCNVlog2FC_thresh)
  idx_padj        = which(!is.na(gene_list$padj       ) &     gene_list$padj         < padj_thresh)  
  idx = intersect(intersect(idx_padj, idx_z_score), idx_noCNVlog2FC)
  candidate_list = gene_list[idx, ]
  }
  return(candidate_list)
}


# select_candidates = function(gene_list, padj_thresh, z_score_thresh, noCNVlog2FC_thresh){
#   idx_padj        = !is.na(gene_list$padj       ) &     gene_list$padj         < padj_thresh       
#   idx_z_score     = !is.na(gene_list$z_score    ) & abs(gene_list$z_score    ) < z_score_thresh    
#   idx_noCNVlog2FC = !is.na(gene_list$noCNVlog2FC) & abs(gene_list$noCNVlog2FC) > noCNVlog2FC_thresh
#   idx = intersect(intersect(idx_padj, idx_z_score), idx_noCNVlog2FC)
#   candidate_list = gene_list[idx, ]
#   return(candidate_list)
# }

  
