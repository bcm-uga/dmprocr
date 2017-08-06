### perform linear regression on data_trscr and data_cnv  
### (be careful to use normalized count data as input)
dmprocr_lin_reg = RNAseq_cnv_reg(data_ntrscr = dmprocr_result$data_ntrscr, 
                                     data_cnv = dmprocr_study$data_cnv,
                                     exp_grp = dmprocr_study$exp_grp, 
                                     gene_list = dmprocr_result$result_list, 
                                     contrast = c("sample","case","ctrl"))

### you can change the filter on CNV distribution 
### and look only at a selected subset of indivuals (here females)
filter_indiv = rownames(dmprocr_study$exp_grp[dmprocr_study$exp_grp$gender == "F", ])
dmprocr_lin_reg_f =  RNAseq_cnv_reg(data_ntrscr = dmprocr_result$data_ntrscr, 
                                 data_cnv = dmprocr_study$data_cnv,
                                 exp_grp = dmprocr_study$exp_grp, 
                                 gene_list = dmprocr_result$result_list,
                                 filter_indiv = filter_indiv,
                                 contrast = c("sample","case","ctrl"),
                                 cnv_filter = c(0.05, 0.95))

