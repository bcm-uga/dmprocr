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
#' @param contrast A vector containing the constrast to be used to estimate the logarithmic fols change
#' @param fitType A DESeq fitting paramater, by default set to "parametric"
#' 
#' @return A \code{gene_list} table including log2FoldChange and adjusted p-value (padj) computed by DESeq2 and a \code{data_ntrscr} matrix of normalized counts.
#' 
#' @example examples/example-generate_fakestudy.R
#' @example examples/example-RNAseq_diffAnalysis.R
#' 
#' @importFrom DESeq2 "DESeqDataSetFromMatrix"
#' @importFrom DESeq2 "DESeq"
#' @importFrom DESeq2 "results"
#' @importFrom DESeq2 "counts"
#' @importFrom stats "as.formula"
#'
#' @export

RNAseq_diffAnalysis <- function(data_trscr, exp_grp, gene_list, filter_indiv = "no_filter", alpha = 0.05, contrast=c("sample","01","11"), fitType="parametric") {
  
  if (filter_indiv[1] == "no_filter") {
    print("all individuals will be used in the differential analysis")
    filter_indiv = colnames(data_trscr)
    }
  if (length(which(filter_indiv %in% colnames(data_trscr))) != length(filter_indiv) ) {
    stop("ERROR, filter_indiv does not fit to indivuals present in the data matrix")
    }
  
  #Generate the dds for DESeq2 (we start from a count matrix)
  countData <- data_trscr[gene_list$gene_id,  filter_indiv] #data matrix
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

  #Perform the differential analysis 
  dds <- DESeq2::DESeq(dds, fitType= fitType)
 
  #Extract the results of differential analysis
  res <- DESeq2::results(dds,alpha = alpha ,  contrast=contrast)
  result_list = cbind(gene_list[,1:7], res[gene_list$gene_id, c("log2FoldChange","padj")])
  
  #Generate the normalized data matrix
  data_ntrscr = DESeq2::counts(dds, normalized=TRUE)
  return(list(result_list = result_list, data_ntrscr = data_ntrscr))
  
}

#---------------------------------------------
#' Perform linear regression on expression and CNV data.
#' 
#' Perform linear regression on expression (counts) and CNV data.
#' @param data_ntrscr A \code{data} matrix that contains normalized RNAseq counts (from DESeq2 analysis). Columns correspond to indivuals, row correspond to genes.
#' @param data_cnv A \code{data} matrix that contains CNV data
#' @param exp_grp A \code{exp_grp} dataframe that contains metadatas on \code{data_trscr} individuals.
#' @param gene_list A \code{gene_list} bedfile containing the genes for which the linear regression will be perform.
#' @param filter_indiv A vector of individual names to be screened for differential expression. Optionnal (set on "no_filter" by default).
#' @param contrast A vector containing the constrast to be used to read metadata
#' @param cnv_filter A vector of two values indicating between which quantiles (for cnv data) the regression should be performed, by default set to c(0.025, 0.975). Should be set to FALSE if no filter is required.
#' 
#' @return A \code{gene_list} table including a z_score value associated with the linear model to the gene_list used in entry.
#'
#' @example examples/example-generate_fakestudy.R
#' @example examples/example-RNAseq_diffAnalysis.R
#' @example examples/example-RNAseq_cnv_reg.R
#' 
#' @importFrom stats "quantile"
#' @importFrom stats "lm"
#'
#' @export

RNAseq_cnv_reg <- function(data_ntrscr, data_cnv, exp_grp, gene_list, filter_indiv = "no_filter", contrast=c("sample","01","11"), cnv_filter = c(0.025, 0.975)){
  
  if (filter_indiv[1] == "no_filter") {
    print("all individuals will be used in the differential analysis")
    filter_indiv = rownames(exp_grp)}
  
  #Select tumor samples for individual that display cnv and trscr data
  cur = exp_grp[exp_grp[[contrast[1]]] == contrast[2], ] # selection of pathogical samples
  tc_indivs = rownames(cur[cur$trscr == 1 & 
                             cur$cnv == 1 &
                             cur$dmprocr_ID %in% filter_indiv, ])
  
  #Generate dataframe for linear regression with a z_score corresponding to the standardized beta coefficient 
  
  if (length(cnv_filter) == 1) {
    print("patient are not filter upon cnv values")
  } else {
    print("patient are filter upon cnv values")
  }
  
  foo = apply (t(rownames(gene_list)), 2, function (gene) {
    
    i = which(rownames(gene_list) == gene)
    if (i %% 100 == 0){print (paste(gene,i, "/", length(rownames(gene_list)), sep=" "))}
    
    #Generate a temporary dataframe with counts values (exp) and cnv value for each patient
    tmp_data = apply (t(tc_indivs), 2, function (patient) {
      exp = data_ntrscr[gene, patient]
      cnv = data_cnv[gene, patient]
      return(list(patient = patient, exp=exp, cnv = cnv))
    })
    cur_data = (do.call(rbind, tmp_data))
    cur_data = data.frame(lapply(data.frame(cur_data, stringsAsFactors=FALSE), unlist), stringsAsFactors=FALSE)
    
    #Check if values are different of )
    if (sum(cur_data$cnv) == 0 ) {
      z_score = NA
    } else {
      
      #If no filter is required, perform the regression on all data
      if (length(cnv_filter) == 1) {
        mod = lm(exp ~ cnv, data=cur_data)
        
      #If filter required, perform the regression on reduced dataframe
      } else {
        red_cur_data = cur_data[cur_data$cnv < quantile(cur_data$cnv, cnv_filter[2]) &
                                  cur_data$cnv > quantile(cur_data$cnv, cnv_filter[1]), ]
        mod = lm(exp ~ cnv, data=red_cur_data)
      }
      
      #Extract the z_score from the model
      z_score = summary(mod)$coefficients[6]
    }
    
    gene_id = gene
    return(list(gene_id = gene_id, z_score = z_score))
   })
  data_reg = do.call(rbind, foo)   
  data_reg = data.frame(lapply(data.frame(data_reg, stringsAsFactors=FALSE), unlist), stringsAsFactors=FALSE)
  
  #merging with existing gene_list
  if (length(which(data_reg$gene_id %in% gene_list$gene_id)) != length(data_reg$gene_id) ) {
    stop("ERROR, some genes are missing in the linear regression z_score")
  }
    
    result_list = merge(gene_list, data_reg, by='gene_id')
    rownames(result_list) = result_list$gene_id
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
#' @param contrast A vector containing the constrast to be used to estimate the logarithmic fols change
#' @param no_cnv_filter A vector of two values indicating the CNV values threshold for no CNV, by default set to c(-0.2, 0.2). 
#' 
#' @example examples/example-generate_fakestudy.R
#' @example examples/example-RNAseq_diffAnalysis.R
#' @example examples/example-RNAseq_cnv_reg.R
#' @example examples/example-noCNV_diffAnalysis.R
#' 
#' @return A \code{gene_list} table including a noCNVlog2FC value associated with differential expression between sample without CNV.
#'
#' @export

noCNV_diffAnalysis = function(data_ntrscr, data_cnv, exp_grp, gene_list, filter_indiv = "no_filter", contrast=c("sample","01","11"), no_cnv_filter = c(-0.2, 0.2)){
  foo = apply (t(rownames(gene_list)), 2, function (gene) {
    i = which(rownames(gene_list) == gene)
    if (i %% 500 == 0){print (paste(gene,i, "/", length(rownames(gene_list)), sep=" "))}
    
    #identify samples with no CNVs
    if (mean(data_cnv[gene,]) == 0) {
      noCNVlog2FC = NA
    } else if (mean(data_cnv[gene, ]) != 0){
      idx_samples = colnames(data_cnv)[data_cnv[gene, ] > no_cnv_filter[1] 
                                       & data_cnv[gene, ] < no_cnv_filter[2]]
      cur = exp_grp[exp_grp$dmprocr_ID %in% idx_samples & exp_grp$trscr == 1, ]
      idx_p = rownames(cur)[cur[[contrast[1]]] == contrast[2]] #look for pathological samples
      idx_c = rownames(cur)[cur[[contrast[1]]] == contrast[3]] #look for ctrl samples
      ntrscr_p = mean(as.numeric(data_ntrscr[gene,idx_p]))
      ntrscr_c = mean(as.numeric(data_ntrscr[gene,idx_c]))
      noCNVlog2FC = log2(ntrscr_p/ntrscr_c)
    }
    return(list(gene_id = gene, noCNVlog2FC = noCNVlog2FC ))
  })
  noCNV_data = do.call(rbind, foo)   
  noCNV_data = data.frame(lapply(data.frame(noCNV_data, stringsAsFactors=FALSE), unlist), stringsAsFactors=FALSE)
  result_list = merge(gene_list, noCNV_data, by='gene_id')
  rownames(result_list) = result_list$gene_id
  return(result_list)
}

#---------------------------------------------
#' Generate a candidate table
#' 
#' @param gene_list A \code{gene_list} bedfile containing the log2FC, zscore and noCNVlog2FC values.
#' @param t_padj A threshold p-value under which expression is significantly different by DESeq analysis
#' @param t_z_score A threshold zscore value, candidate genes should display a zcore value such as -t_zscore<zscore<t_zscore
#' @param t_noCNVlog2FC A threshold noCNVlog2FC, candidate genes should display a noCNVlog2FC value such as noCNVlog2FC > t_noCNVlog2FC or noCNVlog2FC < t_noCNVlog2FC
#' 
#' @example examples/example-generate_fakestudy.R
#' @example examples/example-RNAseq_diffAnalysis.R
#' @example examples/example-RNAseq_cnv_reg.R
#' @example examples/example-noCNV_diffAnalysis.R
#' @example examples/example-select_candidates.R
#' 
#' @return A \code{gene_list} table including of candidate genes
#' @export


select_candidates = function(gene_list, t_padj, t_z_score, t_noCNVlog2FC){
  candidate_list = gene_list[!is.na(gene_list$z_score) & !is.na(gene_list$padj) & !is.na(gene_list$noCNVlog2FC), ]
  candidate_list = candidate_list[candidate_list$z_score > (-t_z_score) & 
                                    candidate_list$z_score < t_z_score &
                                    candidate_list$padj < t_padj &
                           (candidate_list$noCNVlog2FC > t_noCNVlog2FC | candidate_list$noCNVlog2FC < (-t_noCNVlog2FC) )  , ]
  return(candidate_list)
}

  
