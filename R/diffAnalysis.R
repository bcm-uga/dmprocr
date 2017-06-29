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
#' 
#' @return A \code{gene_list} table including log2FoldChange and adjusted p-value (padj) computed by DESeq2 and a \code{data_ntrscr} matrix of normalized counts.
#' 
#' @importFrom DESeq2 "DESeqDataSetFromMatrix"
#' @importFrom DESeq2 "DESeq"
#' @importFrom DESeq2 "results"
#' @importFrom DESeq2 "counts"
#'
#' @export

RNAseq_diffAnalysis <- function(data_trscr, exp_grp, gene_list, filter_indiv = "no_filter", alpha = 0.05, contrast=c("sample","01","11")) {
  
  if (filter_indiv[1] == "no_filter") {
    print("all individuals will be used in the differential analysis")
    filter_indiv = colnames(data_trscr)}
  if (length(which(filter_indiv %in% colnames(countData))) != length(filter_indiv) ) {
    stop("ERROR, filter_indiv does not fit to indivuals present in the data matrix")
    }
  
  #Generate the dds for DESeq2 (we start from a count matrix)
  countData <- data_trscr[gene_list$gene_id,  filter_indiv] #data matrix
  colData <- exp_grp[filter_indiv, c("dmprocr_ID", "sample")]  #sample type 
  
  #Check wether columns ids of the count matrix corresponds to rows ids of the column data 
  if (length(which(rownames(colData) == colnames(countData))) != ncol(countData) ) {
    stop("ERROR, exp_grp does not correspond to data matrix")
    }
  
  #Construct a DESeqDataSet
  colData$sample <- as.factor(colData$sample)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = countData,
                                colData = colData,
                                design = ~ sample)
  
  #Perform the differential analysis 
  dds <- DESeq2::DESeq(dds)
 
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
#' @param gene_list A \code{gene_list} bedfile containing the genes to screen for differential expression.
#' @param filter_indiv A vector of individual names to be screened for differential expression. Optionnal (set on "no_filter" by default).
#' @param tumor A factor indicating the label of pathological samples in the sample column, set to "01" by default
#' @param cnv_filter A vector of two values indicating between which quantiles (for cnv data) the regression should be performed, by default set to c(0.025, 0.975). Should be set to FALSE if no filter is required.
#' 
#' @return A \code{gene_list} table including a z_score value associated with the linear model to the gene_list used in entry.
#'
#' @importFrom stats "quantile"
#' @importFrom stats "lm"
#'
#' @export

RNAseq_cnv_reg <- function(data_ntrscr, data_cnv, exp_grp, gene_list, filter_indiv = "no_filter", tumor = "01", cnv_filter = c(0.025, 0.975)){
  
  if (filter_indiv[1] == "no_filter") {
    print("all individuals will be used in the differential analysis")
    filter_indiv = rownames(exp_grp)}
  
  #Select tumor samples for individual that display cnv and trscr data
  tc_indivs = rownames(exp_grp[ exp_grp$sample == tumor &
                                  exp_grp$trscr == 1 & 
                                  exp_grp$cnv == 1 &
                                  exp_grp$dmprocr_ID %in% filter_indiv, ])
  
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
    return(result_list)
}
