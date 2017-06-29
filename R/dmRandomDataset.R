#'dmRandomDataset
#' 
#'Generate random meth_study to run the functions examples.
#'
#'@example examples/example-dmRandomDataset.R
#'
#'
#'@export
dmRandomDataset = function() {
  # genes
  txStart = unique(round(runif(100)*1000000))
  genes       <- data.frame(chr="chr1", txStart=txStart, txStop=txStart + round(abs(rnorm(length(txStart)) * 10000)), geneSymbol="gene", refGen="simu", 
                            strand=ifelse(runif(length(txStart))>0.5, "+", "-"), stringsAsFactors=FALSE)
  genes = genes[order(genes$chr, genes$txStart, genes$txStop),]
  genes$geneSymbol = paste0("gene_", 1:length(txStart))
  head(genes)
  dim(genes)
  
  # exp_grp
  exp_grp = data.frame(ref=rep(c("case", "ctrl"), 25), 
                       patientId=rep(1:25, each=2), 
                       row.names=paste0("sample_", 1:50), stringsAsFactors=FALSE)
  head(exp_grp)
  dim(exp_grp)
  
  # platform
  probe_pos = apply(genes, 1, function(gene){
    strand = gene[[6]]
    chr = gene[[1]]
    if (strand=="+") {
      tss = as.integer(gene[[2]])
    } else {
      tss = as.integer(gene[[3]])
    }
    unique(round(rnorm(round(runif(1, 10,30)), tss, round(runif(1, 300, 1000)))))
  })
  probe_pos = sort(unique(unlist(probe_pos)))
  platform = data.frame(chr="chr1", pos=probe_pos, row.names=paste0("probe_", 1:length(probe_pos)), stringsAsFactors=FALSE)
  head(platform)
  dim(platform)
  
  # data
  data = matrix(runif(nrow(platform) * nrow(exp_grp),0,1), nrow(platform), nrow(exp_grp))
  rownames(data) = rownames(platform)
  colnames(data) = rownames(exp_grp)
  
  # study
  meth_study = list(data=data, exp_grp=exp_grp, platform=platform)
  return(list(meth_study=meth_study, genes=genes))
}