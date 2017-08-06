### perform DE analysis on individuals with no cnv variation 
### (be careful to use normalized count data as input)
dmprocr_lin_reg_no_cnv = noCNV_diffAnalysis(data_ntrscr = dmprocr_result$data_ntrscr, 
                                     data_cnv = dmprocr_study$data_cnv,
                                     exp_grp = dmprocr_study$exp_grp, 
                                     gene_list = dmprocr_lin_reg, 
                                     contrast = c("sample","case","ctrl"))

### you can change the filter on CNV 
### (considere no cnv variation of cnv values between -0.3 and 0.3) 
### and look only at a selected subset of indivuals (here females)
filter_indiv = rownames(dmprocr_study$exp_grp[dmprocr_study$exp_grp$gender == "F", ])
dmprocr_lin_reg_no_cnv_f =  noCNV_diffAnalysis(data_ntrscr = dmprocr_result$data_ntrscr, 
                                 data_cnv = dmprocr_study$data_cnv,
                                 exp_grp = dmprocr_study$exp_grp, 
                                 gene_list =  dmprocr_lin_reg_f,
                                 filter_indiv = filter_indiv,
                                 contrast = c("sample","case","ctrl"),
                                 no_cnv_filter = c(-0.3, 0.3))

