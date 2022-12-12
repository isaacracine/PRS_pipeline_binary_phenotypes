cov_file <- snakemake@input[[1]]
eigenvec_file <- snakemake@input[[2]]
target_name <- snakemake@params[[1]]

setwd('/staging/leuven/stg_00092/IBP_PRSproject/1_Intestinal_Bowel_Disease/data')

cov <- read.table(cov_file, header=T)
pcs <- read.table(eigenvec_file, header=F)
colnames(pcs) <- c("FID","IID", paste0("PC",1:6))
covariate <- merge(cov, pcs, by=c("FID", "IID"))
covariate_name <- paste(target_name, ".covariate", sep = "")

write.table(cov,snakemake@output[[1]], quote=F, row.names=F)
q()
