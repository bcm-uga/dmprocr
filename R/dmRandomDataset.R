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
  txStart = unique(round(stats::runif(100)*1000000))
  genes       <- data.frame(chr="chr1", txStart=txStart, txStop=txStart + round(abs(stats:: rnorm(length(txStart)) * 10000)), geneSymbol="gene", refGen="simu", 
                            strand=ifelse(stats::runif(length(txStart))>0.5, "+", "-"), stringsAsFactors=FALSE)
  genes = genes[order(genes$chr, genes$txStart, genes$txStop),]
  genes$geneSymbol = paste0("gene_", 1:length(txStart))
  utils::head(genes)
  dim(genes)
  
  # exp_grp
  exp_grp = data.frame(ref=rep(c("case", "ctrl"), 25), 
                       patientId=rep(1:25, each=2), 
                       row.names=paste0("sample_", 1:50), stringsAsFactors=FALSE)
  utils::head(exp_grp)
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
    unique(round(stats:: rnorm(round(stats::runif(1, 10,30)), tss, round(stats::runif(1, 300, 1000)))))
  })
  probe_pos = sort(unique(unlist(probe_pos)))
  platform = data.frame(chr="chr1", pos=probe_pos, row.names=paste0("probe_", 1:length(probe_pos)), stringsAsFactors=FALSE)
  utils::head(platform)
  dim(platform)
  
  # data
  data = matrix(stats::runif(nrow(platform) * nrow(exp_grp),0,1), nrow(platform), nrow(exp_grp))
  rownames(data) = rownames(platform)
  colnames(data) = rownames(exp_grp)
  
  # study
  meth_study = list(data=data, exp_grp=exp_grp, platform=platform)
  return(list(meth_study=meth_study, genes=genes))
}

#---------------------------------------------
#'Generate simulated dataset for dmprocr package.
#'
#'Generate transcriptome (HTseq counts), cnv (copy number segment) and methylation (beta values) for 500 genes, 9810 methylation probes, in 50 samples.
#'
#'@example examples/example-generate_fakestudy.R
#'
#'
#'@export
generate_fakestudy = function() {
  # gene_list
  start = unique(round(stats::runif(500)*1000000))
  gene_list = data.frame(chr="chr1",
                         start =start,
                         stop = start + round(abs(stats:: rnorm(length(start)) * 10000)),
                         gene_id="gene",
                         score = "NA",
                         strand = ifelse(stats::runif(length(start))>0.5, "+", "-"),
                         ref_genome="simu",
                         stringsAsFactors=FALSE)
  gene_list$gene_id = paste0("gene_", 1:length(start))
  rownames(gene_list) = gene_list$gene_id
  utils::head(gene_list)
  dim(gene_list)
  
  
  # exp_grp
  sample = rep(c("case", "ctrl"), 25)
  gender = c(rep(c("F", "F", "M", "M"), 12), "F", "F")
  patient_ID = rep(1:25, each=2)
  dmprocr_ID = paste(patient_ID, sample, sep = "_")
  trscr = rep(1, 50)
  cnv = rep(1, 50)
  meth = rep(1, 50)
  exp_grp = data.frame(dmprocr_ID, patient_ID, sample, trscr, cnv, meth, gender,
                       row.names=dmprocr_ID)
  utils::head(exp_grp)
  dim(exp_grp)
  
  # platform
  probe_pos = apply(gene_list, 1, function(gene){
    strand = gene[[6]]
    chr = gene[[1]]
    if (strand=="+") {
      tss = as.integer(gene[[2]])
    } else {
      tss = as.integer(gene[[3]])
    }
    unique(round(stats:: rnorm(round(stats::runif(1, 10,30)), tss, round(stats::runif(1, 300, 1000)))))
  })
  probe_pos = sort(unique(unlist(probe_pos)))
  pf_meth = data.frame(chr="chr1", pos=probe_pos, row.names=paste0("probe_", 1:length(probe_pos)), stringsAsFactors=FALSE)
  utils::head(pf_meth)
  dim(pf_meth)
  
  # data_trscr
  data_trscr = matrix(stats::rnbinom(n = nrow(gene_list) * nrow(exp_grp), size=100,  mu=800), nrow(gene_list), nrow(exp_grp))
  rownames(data_trscr) = rownames(gene_list)
  colnames(data_trscr) = rownames(exp_grp)
  for (i in 1:dim(data_trscr)[1]){
    if (i %% 25 == 1) { data_trscr[i, seq(1,50,2)] = round(data_trscr[i,seq(1,50,2)]*sample(c(1/(3:10), 3:10),1))} #generates deregulated genes in case samples
  }
  utils::head(data_trscr)
  dim(data_trscr)
  
  # data_cnv
  data_cnv = matrix(stats::rnorm(n = nrow(gene_list) * nrow(exp_grp), mean=0,  sd=0.4), nrow(gene_list), nrow(exp_grp))
  rownames(data_cnv) = rownames(gene_list)
  colnames(data_cnv) = rownames(exp_grp)
  utils::head(data_cnv)
  dim(data_cnv)
  
  # data_meth
  data_meth = matrix(stats::runif(nrow(pf_meth) * nrow(exp_grp),0,1), nrow(pf_meth), nrow(exp_grp))
  rownames(data_meth) = rownames(pf_meth)
  colnames(data_meth) = rownames(exp_grp)
  utils::head(data_meth)
  dim(data_meth)
  
  # study
  return(list(data_trscr = data_trscr,
              data_cnv = data_cnv,
              data_meth = data_meth,
              gene_list = gene_list,
              pf_meth = pf_meth,
              exp_grp = exp_grp))
  
}
