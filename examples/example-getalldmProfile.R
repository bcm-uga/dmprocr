### get dmProfiles for multiple genes
genes        <- randomDataset$genes 
alldmProfile <- getalldmProfile(diff_meth_study,
                                bedfile = genes, slide = 500,
                                pf_chr_colname = "chr",
                                pf_pos_colname = "pos")

