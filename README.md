# <ins>PRS pipeline for binary phenotypes</ins> 
This repository will host the pipeline and necessary files to calculate and compare PRS for multiple tools  (implemented via snakemake)

With this pipeline, we calculate polygenic risk scores (PRS) using four different PRS-calculating tools (PLINK, PRSice, lassosum and LDpred). Additively, we compare the calculated PRS of each method to eachother. 

This project aims to compare four different tools that calculate Polygenic Risk Scores (PRS). These tools are PLINK, PRSice, Lassosum and LDpred. Currently, this workflow only allows for testing of **binary** phenotypes, such as cases and controls.

The heart of this project is a Snakemake script, therefore allowing for extension with more tools. Additionally, the project is supported for users that have access to a remote cluster. The scripts require the submission of a SLURM job, however if your cluster uses another job scheduler this can easily be adapted as well. The aim of using this pipeline is not to take the PRS directly. Rather it should aid the user in selecting the best tool to use for calculating PRS on their specific dataset. The user should rerun their own analysis on such tool with their data to allow for more fine tuning. 

_If you are working in a shared remote directory, remember to give access to the other members after creating any file or directory._

## <ins>PREVIOUS REQUIREMENTS:</ins>

* Conda
* The GWAS base dataset should contain uppercase alleles.
* We require that you have the effect size in beta metric rather than odds ratio. To obtain the beta estimate from an odds ratio we recommend using R. 

```
data <- read.table(<GWAS_summary_statisitcs_file>, header=T)
data$beta <- log(dat$OR)
write.table(data, "<name_of_new_GWAS_summary_statisitcs_file>", quote=F, row.names=F)
```

## <ins>DOWNLOAD THIS REPOSITORY:</ins>

In order to run this tool you have to download this GitHub repository. Then you can modify teh corresponding variables (explained later in this README.md file) and run it locally. In case you wanted to push it to your GitHub repositories, a fork will be automatically created, and if you wanted to propose any change you can submit a pull request.
```
git clone [https://github.com/jcasadogp/IBP_PRS_2022.git](https://github.com/isaacracine/PRS_pipeline_binary_phenotypes.git)
```

It includes the following folders:
* Conda environments (`PRS_pipeline_binary_phenotypes/conda_env/`): several .yml files that will be used by the different scripts to activate the contained conda environments. These .yml files are what allows for the running of snakefile, lassosum, LDpred and generation of performance metric plots.
* Plink (`PRS_pipeline_binary_phenotypes/plink/`): installation and executable files for PLINK.
* Prsice (`PRS_pipeline_binary_phenotypes/prsice/`): installation, executable and R script for PRSice.
* Scripts (`PRS_pipeline_binary_phenotypes/scripts/`): it contains three types of files: R scripts, Snakefiles and the SLURM jobs' files that must be run by the user.


## <ins>OVERVIEW: </ins>

During development we stumbled across a problem implementing LDpred due to its drastically large memory requirements. Therefore, we developed two pipelines within this repository, one with and without LDpred. If the user would like to use the pipeline with LDpred they must also have a either a **small dataset** or access to the **unimputed data**. If either of these expectations are not met then the user will likely not have sufficient resources. 

Therefore, the `scripts`directory contains the following Snakefiles and .pbs files:
* `LD_pred_all_job.pbs`: SLURM job requierd for running ```LD_pred_all_snakefile```
* `LD_pred_all_snakefile`: the snakemake (pipeline) which runs also LDpred
* `imp_only_job.pbs`: SLURM job required for running ```imp_only_snakefile```
* `imp_only_snakefile`: the snakemake (pipeline) which does not run LDpred


## <ins>INPUT DATA:</ins>

Data needed:

* GWAS summary statistics
* 1000 genome file to use as reference, which should be a subset of the cohort that is of the same ancestral background as your target cohort
* Target files: PLINK binary format files (.fam, .bed, .bim)
* Phenotype file: File containing the phenotype file with FID, IID, and the phenotype must be encoded as **1 for controls and 2 for cases**.
* Covariates and eigenvectos files (eigenvalues are not really necessary). To generate the eigenvector file you can must perform pruning first. Use the following commands:
```
plink \
    --bfile <target_data_prefix> \
    --indep-pairwise 200 50 0.25 \ #can change these values (look at PLINK documentation for further explaination)
    --out <target_data_prefix>
plink \
    --bfile <target_data_prefix> \
    --extract <target_data_prefix>.prune.in \
    --pca 6 \ #we only consider the first 6 PCs in the pipeline, retaining fewer will cause errors
    --out <target_data_prefix>
```

All this data will be placed in the `data/` directory. To create it run:

```
cd PRS_pipeline_binary_phenotypes/
mkdir data/
cp <target_files> data/
```

## <ins>MAKE A SNAKEMAKE CONDA ENVIRONMENT:</ins>

We require the user to create a conda environment for running snakemake. This is to prevent issues with different versions of snakemake.

```
cd PRS_pipeline_binary_phenotypes/conda_env/
conda env create -f snakemake_PRS.yml
```

Important: note down the directory in which the environment was created, you will need it to adapt the slurm job script. This is often within your miniconda directory. 


## <ins>OUTPUT DIRECTORY:</ins>

The output data will be placed in a directory called ```output_data/``` that contains a different directory per tool. Also, inside each of the tools' directories, we will separate the data into ```target_data/``` and ```external_data/```, depending on if it comes from the tool with the 1000genome data (```external_data/```) or from the target data (```target_data/```). To create all these directories run:

```
cd IBP_PRS_2022
mkdir output_data/
cd output_data/
mkdir -p 001_plink/
mkdir -p 001_plink/target_data/
mkdir -p 001_plink/external_data/
mkdir -p 002_prsice/
mkdir -p 002_prsice/target_data/
mkdir -p 002_prsice/external_data/
mkdir -p 003_lassosum/
mkdir -p 003_lassosum/target_data/
mkdir -p 003_lassosum/external_data/
mkdir -p 004_LDpred/
mkdir -p 004_LDpred/target_data/
mkdir -p 004_LDpred/external_data/
mkdir -p 005_comparison/
```

## <ins>ADAPTING THE SNAKEMAKE FILE:</ins>
Snakemake is a wonderful workflow engine which allows for easy adaptability and extension. At the top of the script you will see several global variables, all of which you will need to fill in. They include the desired reference genome for a list of options, different p-value thresholds to test in PLINK and PRSice, different shrinkage values to test in lassosum and the names of specific columns from the GWAS summary statistic file. Here is an example:

imp_only_snakefile:
```
project_dir = <parent-directory> + "/PRS_pipeline_binary_phenotypes/" #Change <parent-directory> by the directory in which you cloned the repository

# === Prefix of the files  ===
TARGET_files = "IBD_GSA_fin_maf" #prefix of target file
EXTERNAL_file = "1000G_EUR_fin_maf" #prefix of external (1000 Genomes) file
GWAS_file = "GWAS_summary_stats.txt" #summary statistics file
PHENO_file = "final_phenotypes.txt" #phenotype file
# === Parameters for tools ===
plink_pval_thresholds = "5e-8,1e-5,0.01,0.05,0.1,0.5" #must be in a comma separated format with values from 0 to 1, NO SPACES
lasso_thresholding_values = "0.2,0.5,0.9,1" #must be in a comma separated format with values from 0 to 1, NO SPACES
genome_build = "EUR.hg38" #options: "EUR.hg19", "AFR.hg19", "ASN.hg19", "EUR.hg38", "AFR.hg38", "ASN.hg38"
noneffect_allele = "Allele1" 
effect_allele = "Allele2" #the allele which the beta is associated with 
chromosome_col = "CHR"
pos_col = "pos_hg38" #the column of the basepair positions of the variants
beta_col = "Effect" #the column of betas, effect size, of the variants
p_value_col = "P.value" #the column of associated p-values of the variants 
snp_col = "SNP_hg38" #the column of the SNP identifiers of the variants
maf_value = "0.01" #can be a MAF threshold or the name of the column containing MAF estimates
info_value = "0.8"
binary_pheno = "T" #indicate if the phenotype is a binary trait or not. Values: "T" or "F"
```

LD_pred_all_snakefile: in addition to the above specified variables you will also have to adapt the following variables for running LDpred:
```
TARGET_ldpred_files = "IBD_GSA_fin_maf" #prefix of the LDpred target file
EXTERNAL_ldpred_files = "1000G_EUR_fin_maf" #prefix of external 1000 Genomes file
LDpred_num_pvales = 1 #must be a positive integer, but we recommend no more than 4 for run time 
```


The user should provide the file names and column names as above. The thresholding and shrinkage parameters can be adapted to allow for more niche testing.

The more thresholds, shrinkage parameters and p-values for LDpred that are tested the longer the pipeline will need. It's estimated that the ```imp_only_snakefile``` script will take around 6 hours, while the ```LD_pred_all_snakefile``` will likely take more than 24 hours to run.

## <ins>ADAPTING THE SLURM JOB SCRIPT:</ins>

This pipeline was developed using the Vlaams Supercomputer Centrum and utilizes a SLURM job scheduler. The user should look into details about their own institution's cluster and scheduler. Begin your BASH script by setting the directory to the global snakemake file downloaded from this repository. Further, makesure you have created an envinroment from the `snakemake.yml` we provided as this will prevent version and dependcy issues.

Important: do not place the pbs and the snakefile on separated directories. If you do so, add a ```cd <directory_to_snakemake>``` line in the pbs file so it finds the Snakefile to be run.

<ins>imp_only_job.pbs: </ins>
We recomend the user allote at least 10 hours for the pipeline to run as the shrinkage based methods are rather computationally intesive. This will avoid running into walltime issues. Further, these shrinkage tools require large amounts of data so do not be shy with the amount of data requested.


```
#!/usr/bin/env bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=10:00:00
#PBS -l pmem=40gb
#PBS -A <account>

#cd <directory_to_snakemake>
module purge
eval "$(conda shell.bash hook)"
conda activate snakemake_PRS
snakemake -s imp_only_snakefile --cores 4 --use-conda
```

<ins>LD_pred_all_job.pbs: </ins>
To run the snakemake considering LDpred we strongly beleive access to a partition on your cluster that allows for very large memory allocation is needed, hence why we specified `#PBS -l pmem=720gb` and `#PBS -l partition = bigmem`. Therefore, it is recommended to be aggressive with the amount of memory and the allocated wall time. Besides the header take care to pay attention to the name of the snakefile has also changed.
```
#!/usr/bin/env bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=30:00:00
#PBS -l pmem=720gb
#PBS -l partition = bigmem
#PBS -A <account>

#cd <directory_to_snakemake>
module purge
eval "$(conda shell.bash hook)"
conda activate snakemake_PRS
snakemake -s LD_pred_all_snakefile --cores 4 --use-conda
```

Results of the pipeline will be written to log files. Check the error file to see if all the steps in the pipeline completed. If you have an error please first familarize yourself with Snakemake (see sources). 

Once this is set up **CONGRADULATIONS**! You can now submit this job to your cluster's scheduler and patiently wait for your results. 

## <ins>OUTPUT FILES:</ins>
Several intermediate files, files containing statistics and predicted PRS will be generated. Once again the user should focus their attention to the `005_comparison/` directory. This is because this directory will contain the useful plots comparing the performance of the tools. It will also contain a .txt file consisiting of the best parameter values selected for each tool, `best_thr_all.txt`. Once again, the PRS generated from these tools are **not** meant to be used directly for further analyses or publIshing, rather to help the user to select which tool will have the best performance on their dataset.

## <ins>OUTPUT PLOTS:</ins>
Several plots will be generated for the prediction of PRS and comparison of the performance metrics across tools.
* A boxplot containing the four tools will be generated comparing the cases and controls
* A ROC curve will be generated to show how well the logistic model built from each tools PRS predictions is at correctly classifying cases and controls
* Bar plots comparing the AUC and R<sup>2</sup> values across the tools, generated using 1000 bootstrap samples

## <ins>WARNINGS:</ins>

### LDpred:
Besides the issues with memory and imputation LDpred also has three different version available to run. We selected for the grid model because we can test specific p-value thresholds. Other implementations can be extended by adaptin the scripts. However, to not the auto model offered with take exceptionally long to time run.

### MISSING PHENOTYPES:
PLINK binary files use -9 to represent missing values, however these will not be considered as missing values in our pipeline given some of the tools used do not recognize -9 as missing. It is important to remove these individuals or phenotypes prior to running the pipeline. 

### ADAPTING YOUR LOGISIC MODELS:
When calculating logistic models for explaining the phenotype from PRS results we do not account for principal components or covariates. However, they are passed into the tools when calculating PRS. Our scripts only consider 6 PCs, therefore calculating more prior to the pipeline execution is okay but calculating fewer is **not**. TTherefore, if you wish to include covariates in your logisitc regression you will have to manually adapt the R scripts. One R script, `global_best_plots.R`, decides the best PRS generated across the different thresholds tested. Another R script, `global_performance_plots.R`, produces the plots for comparing the tools. The principal components and covariates can be added by adapting the logisitc regression code as follows:

```
#original command
glm(phenos ~ <prs_column>, data = std_prs_ph.prsice, family = binomial(link = "logit"))
#updated command 
glm(phenos ~ <prs_column> + <pc_1> + <pc_2> + ... + <covariate_1>  + ... + <covariate_n>, data = std_prs_ph.prsice, family = binomial(link = "logit"))
```


# <ins>SOURCES AND FURTHER POINTS OF REFERENCE:</ins> 

If you find yourself stuck please refer to the necessary literature, repositories and websites. These will provide explainations for most of the errors obtainedduring the analysis. Although this is meant to be an automated pipeline that does not mean one can use the pipeline without any thought, so be *cautious*. 

Berisa, T. & Pickrell, J. K. Approximately independent linkage disequilibrium blocks in human populations. Bioinformatics 32, 283-285 (2015).

Choi, S.W., Mak, T.SH. & O’Reilly, P.F. Tutorial: a guide to performing polygenic risk score analyses. Nat Protoc 15, 2759–2772 (2020). 
* User friendly website and guide can be found [here](https://choishingwan.github.io/PRS-Tutorial/)
* Repository can be found [here](https://github.com/choishingwan/PRS-Tutorial)

Mak, Timothy Shin Heng, et al. "Polygenic scores via penalized regression on summary statistics." Genetic epidemiology 41.6 (2017): 469-480. 
* lassosum's repository can be found [here](https://github.com/tshmak/lassosum) repository

 Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, Bender D, Maller J, Sklar P, de Bakker PIW, Daly MJ & Sham PC (2007) 
 PLINK: a toolset for whole-genome association and population-based 
	linkage analysis. American Journal of Human Genetics, 81.
* We utilize PLINK version 1.9
* Link to plinks two websites: [zzz](http://pngu.mgh.harvard.edu/purcell/plink/) and [cog](https://www.cog-genomics.org/plink/1.9/)


Jack Euesden, Cathryn M. Lewis, Paul F. O’Reilly, PRSice: Polygenic Risk Score software, Bioinformatics, Volume 31, Issue 9, 1 May 2015, Pages 1466–1468

Shing Wan Choi, Paul F O'Reilly, PRSice-2: Polygenic Risk Score software for biobank-scale data, GigaScience, Volume 8, Issue 7, July 2019
* A link to the PRSice [repository](https://github.com/choishingwan/PRSice/blob/master/docs/step_by_step.md)

Florian Privé, Julyan Arbel, Bjarni J Vilhjálmsson, LDpred2: better, faster, stronger, Bioinformatics, Volume 36, Issue 22-23, 1 December 2020, Pages 5424–5431
* A link to the [tutorial](https://privefl.github.io/bigsnpr/articles/LDpred2.html)
* A link to the [bigsnpr package](https://www.rdocumentation.org/packages/bigsnpr/versions/1.11.6)


Mölder, F., Jablonski, K.P., Letcher, B., Hall, M.B., Tomkins-Tinch, C.H., Sochat, V., Forster, J., Lee, S., Twardziok, S.O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., Köster, J., 2021. Sustainable data analysis with Snakemake. F1000Res 10, 33. 
* Documentation can be found [here](https://snakemake.readthedocs.io/en/stable/index.html)
