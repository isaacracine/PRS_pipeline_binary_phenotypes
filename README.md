# PRS_pipeline_binary_phenotypes
This repository will host the pipeline and necessary files to calculate and compare PRS for multiple tools  (implemented via snakemake)

With this pipeline, we calculate polygenic risk scores (PRS) using four different PRS-calculating tools (PLINK, PRSice, LDPred and lassoSum). Additively, we compare the calculated PRS of each method to eachother. 


#base dataset (GWAS)

- should contain alleles in uppercase letter. 
- make sure to know which allele is the effect allele. In our case the effect allele is the A2, the non-affected is A1. 



