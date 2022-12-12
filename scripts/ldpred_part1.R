# Load different libraries.
library(remotes)

remotes::install_github("privefl/bigsnpr")
remotes::install_github("privefl/bigstatsr")


library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

library(bigstatsr)

library(magrittr)
library(fmsb)
library(data.table)
library(runonce)

# Read in all the different files.
sumstats <- fread(snakemake@input[[1]])
covariate <- fread(snakemake@input[[2]])
pcs <- fread(snakemake@input[[3]])
phenotype <- fread(snakemake@input[[4]])

snp_readBed(snakemake@input[[5]])
