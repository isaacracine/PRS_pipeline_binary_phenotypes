#this is the script for running a LassoSum analysis on the target
#database
#the install command needs to be run every new interactive session
library(devtools)
install_github("tshmak/lassosum")

library(dplyr)
library(lassosum)
library(data.table)
library(methods)
library(magrittr)
library(parallel)

#invoke 2 threads
#cl <- makeCluster(2)

#read in the files
sum.stat <- snakemake@input[[1]]
bfile <- snakemake@params[[1]]

# Read in and process the covariates
covariate <- fread(snakemake@input[[2]]) 
pcs <- fread(snakemake@input[[3]]) 

#index pcs so only have 8 columns, FID, IID and 6 PCs
pcs <- pcs[,1:8]

#then set the column names 
setnames(pcs, colnames(pcs), c("FID","IID", paste0("PC",1:6)))
#covariate$FID <- as.character(covariate$FID)
#covariate$IID <- as.character(covariate$IID)

# Need as.data.frame here as lassosum doesn't handle data.table 
cov <- merge(covariate, pcs, by = c("FID", "IID"))

#define the LD panel to use
ld.file <- snakemake@params[[3]]

# Read in the target phenotype file
target.pheno <- fread(snakemake@input[[4]])


#our names are messed up, rename
colnames(target.pheno) = c("FID", "IID", "pheno")

#read in the thresholds to test
thresh <- snakemake@params[[2]]

#will need to parse these
thresh.vec <- unlist(strsplit(thresh, ","))
thresh.num <- as.vector(thresh.vec, "numeric")

# Read in the summary statistics
ss <- fread(sum.stat)

########modifying our data format
##need to add a column of N the number of variants
#used in calculating the effect size estimates
#for now I guess just use the number of variants 
ss$N = nrow(ss)

#need to add a column for OR
#need to take the 'effect' column and raise it
#to the exponential 
ss$OR = exp(pull(ss, snakemake@params[[8]]))

#########
# Remove P-value = 0, which causes problem in the transformation
ss <- ss[!snakemake@params[[9]] == 0]

# Transform the P-values into correlation
cor <- p2cor(p = pull(ss, snakemake@params[[9]]),
        n = ss$N, 
        sign = log(ss$OR)
        )

#obtain fam file
fam <- fread(paste0(bfile, ".fam"))
fam[,ID:=do.call(paste, c(.SD, sep=":")),.SDcols=c(1:2)]

# Run the lassosum pipeline
# The cluster parameter is used for multi-threading
# You can ignore that if you do not wish to perform multi-threaded processing
out <- lassosum.pipeline( 
    cor = cor,
    chr = pull(ss, snakemake@params[[6]]),
    pos = pull(ss, snakemake@params[[7]]),
    s = thresh.num,
    A1 = pull(ss, snakemake@params[[5]]),
    A2 = pull(ss, snakemake@params[[4]]),
    ref.bfile = bfile,
    test.bfile = bfile,
    LDblocks = ld.file
    #cluster=cl
)
# we save the validation output
# we use validation rather than pseudovalidation, 
# because we want to allow for covariates 
target.res <- validate(out, pheno = as.data.frame(target.pheno), covar=as.data.frame(cov), plot = FALSE)

#save the best best s, lambda and pgs 
best_s <- target.res$best.s
best_lambda <- target.res$best.lambda
results_pgs <- target.res$results.table

#save the validation table from the results
validation_table <- target.res$validation.table

#save the validation plot
pdf(snakemake@output[[1]]) 
#Model validation refers to the process of confirming that the model actually achieves its intended purpose
plot(target.res)#, ylab = "Validation (correlation)")
# Close the pdf file
dev.off() 
#the best s and lambda will correspond to the highest validation
#value 

#save the result from the validation
write.table(best_s, quote = FALSE, row.names = FALSE, col.names = FALSE, snakemake@output[[2]])
write.table(best_lambda, quote = FALSE, row.names = FALSE, col.names = FALSE, snakemake@output[[3]])
write.table(results_pgs, quote = FALSE, row.names = FALSE, snakemake@output[[4]])
write.table(validation_table, quote = FALSE, row.names = FALSE, snakemake@output[[5]])


###obtain the prs for the best 4 combinations of the 
#thershold parameter, S, and shrinkage parameter, lambda
validation_table$Rank <- rank(-validation_table$value)
best_params <- validation_table[order(validation_table$Rank),]
write.table(best_params, quote = FALSE, row.names = FALSE, snakemake@output[[6]])


l1 <- best_params$lambda[1]
s1 <- best_params$s[1]
l2 <- best_params$lambda[2]
s2 <- best_params$s[2]
l3 <- best_params$lambda[3]
s3 <- best_params$s[3]
l4 <- best_params$lambda[4]
s4 <- best_params$s[4]

out1 <- subset(out, s = s1, lambda = l1)
v1 <- validate(out1, test.bfile = bfile, pheno = as.data.frame(target.pheno), covar=as.data.frame(cov), plot = FALSE)

out2 <- subset(out, s = s2, lambda = l2)
v2 <- validate(out2, test.bfile = bfile, pheno = as.data.frame(target.pheno), covar=as.data.frame(cov), plot = FALSE)

out3 <- subset(out, s = s3, lambda = l3)
v3 <- validate(out3, test.bfile = bfile, pheno = as.data.frame(target.pheno), covar=as.data.frame(cov), plot = FALSE)

out4 <- subset(out, s = s4, lambda = l4)
v4 <- validate(out4, test.bfile = bfile, pheno = as.data.frame(target.pheno), covar=as.data.frame(cov), plot = FALSE)

write.table(v1$results.table, quote = FALSE, row.names = FALSE, col.names = TRUE, snakemake@output[[7]])
write.table(v2$results.table, quote = FALSE, row.names = FALSE, col.names = TRUE, snakemake@output[[8]])
write.table(v3$results.table, quote = FALSE, row.names = FALSE, col.names = TRUE, snakemake@output[[9]])
write.table(v4$results.table, quote = FALSE, row.names = FALSE, col.names = TRUE, snakemake@output[[10]])

# Get the maximum R2
r2 <- max(target.res$validation.table$value)^2
write.table(r2, quote = FALSE, row.names = FALSE, col.names = FALSE, snakemake@output[[11]])
