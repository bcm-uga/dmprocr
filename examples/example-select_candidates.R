### Select candidates based on thresholds for padj, z_score and noCNVlog2FC 
candidates = select_candidates(gene_list = dmprocr_lin_reg_no_cnv,
                               padj_thresh = 0.001, 
                               z_score_thresh = 0.3, 
                               noCNVlog2FC_thresh = 2)



