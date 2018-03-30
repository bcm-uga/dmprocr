#' Analyse transcriptome and methylome study 
#'
#' This function analyse transcriptome and methylome study

#' @param gene_symbol the symbol of the gene  
#' @param s_cnv cnv study
#' @param s_meth  methylome study
#' @param s_trscr transcriptome study
#' @param gene_indexed_probes a list of probes, indexed by gene_symbol
#' @export
trscr_meth_analysis = function(gene_symbol, s_cnv, s_meth, s_trscr, gene_indexed_probes) {
  print(gene_symbol)
  
  meth_probe_idx = gene_indexed_probes[[gene_symbol]]
  if (length(meth_probe_idx) <= 1) {
    return(NULL)
  } 
  meth_data = s_meth$data[meth_probe_idx, ]
  meth_data = meth_data[,apply(is.na(meth_data), 2, sum) / nrow(meth_data) < 0.5]
  meth_data = meth_data[apply(is.na(meth_data), 1, sum) / nrow(meth_data) < 0.5,]
  meth_probe_idx = rownames(meth_data)

  # idx_sample = epimedtools::intersect_rec(colnames(s_trscr$data)[order(s_trscr$data[gene_symbol,])], colnames(meth_data), colnames(s_cnv$data)[abs(s_cnv$data[gene_symbol,])<0.2])
  if (!is.null(s_cnv)) {
    idx_sample = intersect(intersect(colnames(s_trscr$data)[order(s_trscr$data[gene_symbol,])], colnames(meth_data)), colnames(s_cnv$data)[abs(s_cnv$data[gene_symbol,])<0.2])   
  } else {
    idx_sample = intersect(colnames(s_trscr$data)[order(s_trscr$data[gene_symbol,])], colnames(meth_data))      
  }

  # idx_sample = epimedtools::intersect_rec(colnames(s_trscr$data)[order(s_trscr$data[gene_symbol,])], colnames(meth_data))

  if (length(idx_sample) <= 1) {
    return(NULL)
  } 
  meth_data = meth_data[, idx_sample]
  trscr_data = s_trscr$data[gene_symbol,idx_sample]

  scores = sapply(meth_probe_idx, function(cpg_probe) {
    # print(cpg_probe)
    idx = !is.na(meth_data[cpg_probe, idx_sample])
    if (sum(idx) != 0) {
      foo = try(entropy::discretize2d(s_trscr$data[gene_symbol,idx_sample[idx]], meth_data[cpg_probe, idx_sample[idx]], 5, 5))
      if (attributes(foo)$class == "try-error") {
        return(NA)          
      }
      mi = suppressWarnings(entropy::mi.plugin(foo))
      return(mi)
    } else {
      return(NA)
    }
  })
    
  ret = list(
    gene_symbol = gene_symbol,
    # tissue_status = s_trscr$exp_grp[idx_sample,]$tissue_status,
    meth_data = meth_data,
    trscr_data = trscr_data,
    study_name = s_trscr$stuffs$name, 
    scores = scores
  )
  

  return(ret)
}


#' Plot transcriptome and methylome co-analyse
#'
#' This plots transcriptome and methylome co-analyse
#' @param results output of trscr_meth_analysis function
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics arrows
#' @importFrom graphics axis
#' @importFrom graphics image
#' @importFrom graphics layout
#' @importFrom graphics text
#' @importFrom stats sd
#' @export
plot_trscr_meth_analysis = function(results) {
  gene_symbol   = results$gene_symbol   
  meth_data     = results$meth_data     
  trscr_data    = results$trscr_data    
  study_name    = results$study_name    
  scores        = results$scores        
  if (!is.null(results$tissue_status)) {
    tissue_status = results$tissue_status
  } else {
    tissue_status = 1    
  }

  layout(matrix(c(1,1, 4, 2, 2, 3, 2, 2, 3), 3), respect=TRUE)
  main = main = paste0(study_name, " - ", gene_symbol)
  plot(trscr_data, 1:length(trscr_data), main=main, xlab="log2(normalized expression)", ylab=paste0(length(trscr_data), " samples"), yaxt="n", col=as.factor(tissue_status))

  colors=c("royalblue", "springgreen", "yellow", "red")
  colors=c("green", "black", "red")
  cols = colorRampPalette(colors)(20)
  breaks = seq(0,1,length.out=length(cols)+1)

  main = main = paste0("Methylome", " - ", gene_symbol)
  image(meth_data, col=cols, breaks=breaks, xaxt="n", yaxt="n", main=main)
  axis(1, (1:nrow(meth_data)-1)/(nrow(meth_data)-1), rownames(meth_data), las=2)

  # layout(matrix(1:27, 3), respect = TRUE)
  plot(scores, ylab="mi score", type="l", xaxt="n", xlab="")
  axis(1, 1:nrow(meth_data), rownames(meth_data), las=2)

  # gene = genes[gene_symbol,]
  # strand = gene[[6]]
  # if (strand == "+") {
  #   tss = gene[[2]]
  #   beg = gene[[2]]
  #   end = gene[[3]]
  # } else {
  #   tss = gene[[3]]
  #   beg = gene[[3]]
  #   end = gene[[2]]
  # }

  # y_base = 0
  # pos_x = pf_orig[meth_probe_idx, pf_pos_colname]


  # plot(0, 0, col = 0, xlim = c(tss - max(up_str, dwn_str),
  #     tss + max(up_str, dwn_str)), ylim = c(-0.5, 1), ylab = "beta",
  #     main = paste0(gene_symbol), xaxt="n", yaxt="n")
  # axis(1, c(tss-max(up_str, dwn_str), tss, tss+max(up_str, dwn_str)), c(paste0("TSS-", max(up_str, dwn_str)/1000, "kb"), "TSS", paste0("TSS+", max(up_str, dwn_str)/1000, "kb")))
  # axis(2, 0:1)
  #
  # x_pb = seq(tss - max(up_str, dwn_str),  tss + max(up_str, dwn_str), length.out=length(meth_probe_idx))
  # arrows(pos_x, -0.1, pos_x, -0.1+0.01, length = 0)
  # arrows(x_pb, -0.2, x_pb, -0.2+0.01, length = 0)
  # arrows(pos_x, -0.1, x_pb, -0.2+0.01, length = 0)
  # text(pos=2, x_pb+(x_pb[2]-x_pb[1])/2, -0.2-0.01, meth_probe_idx, las=2, srt=90, cex=.3)
  #
  #
  # polygon(c(beg, end, end, beg), c(y_base - 0.01, y_base -
  #     0.01, y_base + 0.01, y_base + 0.01), col = 1)
  #
  # arr_len = (end - beg)/3
  # arr_len = 500
  # if (strand == "-") {
  #     arrows(end, y_base, end, y_base + 0.05, length = 0)
  #     arrows(end, y_base + 0.05, end - arr_len, y_base + 0.05,
  #         length = 0.03)
  # } else {
  #     arrows(beg, y_base, beg, y_base + 0.05, length = 0)
  #     arrows(beg, y_base + 0.05, beg + arr_len, y_base + 0.05,
  #         length = 0.03)
  # }
  #
  # m = apply(meth_data, 1, mean)
  # s = apply(meth_data, 1, sd)
  # lines(pos_x, m)
  # lines(pos_x, m-s, lty=2)
  # lines(pos_x, m+s, lty=2)
}

