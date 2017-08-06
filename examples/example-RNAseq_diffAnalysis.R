### perform differential analysis of data_trscr data
### fitType is set to "mean" because the simulated dataset is small (i.e. 500 genes)
dmprocr_result = RNAseq_diffAnalysis(data_trscr = dmprocr_study$data_trscr, 
                                     exp_grp = dmprocr_study$exp_grp, 
                                     gene_list = dmprocr_study$gene_list, 
                                     alpha = 0.05, 
                                     contrast = c("sample","case","ctrl"),
                                     fitType="mean")

### perform differential analysis of data_trscr data only of female patient
### fitType is set to "mean" because the simulated dataset is small (i.e. 500 genes)
filter_indiv = rownames(dmprocr_study$exp_grp[dmprocr_study$exp_grp$gender == "F", ])
dmprocr_result_f = RNAseq_diffAnalysis(data_trscr = dmprocr_study$data_trscr, 
                                     exp_grp = dmprocr_study$exp_grp, 
                                     gene_list = dmprocr_study$gene_list,
                                     filter_indiv = filter_indiv,
                                     alpha = 0.05, 
                                     contrast = c("sample","case","ctrl"),
                                     fitType="mean")

