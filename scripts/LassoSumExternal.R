#this script will calculate the PRS for the external dataset, often
#the 1000 genomes, to use for standardization of scores

#load packages
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

#define the LD panel to use
ld.file <- snakemake@params[[3]]

#read in the thresholds to test
thresh <- snakemake@params[[2]]

#will need to parse these
thresh.vec <- unlist(strsplit(thresh, ","))
thresh.num <- as.vector(thresh.vec, "numeric")

#read in the validation table obtained from
#the target data
valid.tab.ranked <- fread(snakemake@input[[2]], header = TRUE)

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

#we need to select the PRS scores corresponding with the 
#best parameters from the target dataset
l1 <- as.numeric(valid.tab.ranked$lambda[1])
s1 <- as.numeric(valid.tab.ranked$s[1])
l2 <- as.numeric(valid.tab.ranked$lambda[2])
s2 <- as.numeric(valid.tab.ranked$s[2])
l3 <- as.numeric(valid.tab.ranked$lambda[3])
s3 <- as.numeric(valid.tab.ranked$s[3])
l4 <- as.numeric(valid.tab.ranked$lambda[4])
s4 <- as.numeric(valid.tab.ranked$s[4])

out1 <- subset(out, s = s1, lambda = as.character(l1))
v1 <- validate(out1, test.bfile = bfile, plot = FALSE)

out2 <- subset(out, s = s2, lambda = as.character(l2))
v2 <- validate(out2, test.bfile = bfile, plot = FALSE)

out3 <- subset(out, s = s3, lambda = as.character(l3))
v3 <- validate(out3, test.bfile = bfile, plot = FALSE)

out4 <- subset(out, s = s4, lambda = as.character(l4))
v4 <- validate(out4, test.bfile = bfile, plot = FALSE)


#write external data's PRS scores 
write.table(v1$results.table, quote = FALSE, row.names = FALSE, col.names = TRUE, snakemake@output[[1]])
write.table(v2$results.table, quote = FALSE, row.names = FALSE, col.names = TRUE, snakemake@output[[2]])
write.table(v3$results.table, quote = FALSE, row.names = FALSE, col.names = TRUE, snakemake@output[[3]])
write.table(v4$results.table, quote = FALSE, row.names = FALSE, col.names = TRUE, snakemake@output[[4]])