### get the gene probes information
genes = meth_study$genes[2,]
pf_chr_colname = "chr"
pf_pos_colname = "pos"
gene_study <- getMethylInfo(diff_meth_study,
                            bedline = genes,
                            pf_chr_colname = "chr",
                            pf_pos_colname = "pos")
