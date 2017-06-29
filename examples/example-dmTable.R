### get a methylation differential table
meth_study          <- randomDataset$meth_study
tumoral_ref         <- rownames(meth_study$exp_grp)[meth_study$exp_grp$ref=="case"]
control_ref         <- rownames(meth_study$exp_grp)[meth_study$exp_grp$ref=="ctrl"]
diff_meth_study     <- dmTable(meth_study$data,
                              meth_study$platform,
                              meth_study$exp_grp,
                              tumoral_ref,
                              control_ref)
