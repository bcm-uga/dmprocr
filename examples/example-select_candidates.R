### Select candidates based on thresholds for padj, z_score and noCNVlog2FC previously calculated
candidates = select_candidates(gene_list = dmprocr_lin_reg_no_cnv,
                               t_padj = 0.001, 
                               t_z_score = 0.3, 
                               t_noCNVlog2FC = 2)



