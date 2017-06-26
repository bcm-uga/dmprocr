### generate random samples ids
sample_id <- paste0(rep("sample_", 10), 1:10)
probe_id <- paste0(rep("probe_", 10), 1:10)

### generate data
data <- matrix(runif(100,0,1), 10, 10)
colnames(data) <- sample_id
rownames(data) <- probe_id
head(data)
dim(data)

### generate exp_grp
exp_grp <- data.frame(ref=rep(c("case", "ctrl"), 5), 
                      patientId=rep(1:5, each=2), 
                      row.names=sample_id)
head(exp_grp)
dim(exp_grp)

### generate platfrom
platform <- as.data.frame(matrix(rep("x", 50),  10, 5), row.names=probe_id)
head(platform)
dim(platform)

### get indx of ref
tumoral_ref <- rownames(exp_grp[exp_grp$ref == "case", ])
control_ref <- rownames(exp_grp[exp_grp$ref == "ctrl", ])

