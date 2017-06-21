###generate random samples ids
sampleId <- paste0(rep("sample", 10), 1:10)
idProbes <- paste0(rep("probe", 10), 1:10)

###generate data
data <- matrix(runif(100,0,1), 10, 10)
colnames(data) <- sampleId
rownames(data) <- idProbes

###generate exp_group
exp_group <- as.data.frame(matrix(rep("x", 40),  10, 4), row.names = sampleId)
exp_group$ref <- rep(c("01","11"),5)

###generate platfrom
platform  <-as.data.frame(matrix(rep("x", 50),  10, 5), row.names = idProbes)

###get indx of ref
TumoralRef <- rownames(exp_group[exp_group$ref == "01", ])
ControlRef <- rownames(exp_group[exp_group$ref == "11", ])

