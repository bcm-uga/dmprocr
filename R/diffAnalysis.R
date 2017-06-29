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
#' @param filter_indiv A vector of individual names to be screened for differential expression. Optionnal (set on FALSE by default).
#' @param alpha A parameter to indicate the significance cutoff used by \code{DESeq2::results} funtion for optimizing the independent filtering (by default 0.05). If the adjusted p-value cutoff (FDR) will be a value other than 0.1, alpha should be set to that value.
#' @param contrast A vector containing the constrast to be used to estimate the logarithmic fols change
#' 
#' 
#' @return A \code{gene_list} table include log2FoldChange and adjusted p-value (padj) computed by DESeq2 and a \code{data_ntrscr} matrix of normalized counts.
#' 
#' @importFrom DESeq2 "DESeqDataSetFromMatrix"
#' @importFrom DESeq2 "DESeq"
#' @importFrom DESeq2 "results"
#' @importFrom DESeq2 "counts"
#'
#'
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

