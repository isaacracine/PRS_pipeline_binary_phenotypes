## TO DO

print("Entering the R script: global_plot_snakemake.R")

library(dplyr)
library(ggplot2)
library(data.table)
library(ROCR)
library(pROC)
library(pscl)
library(stringr)

# =========================================================================================================================================================================================

# # GET THEM FROM SNAKEMAKE

project_dir = snakemake@params[[1]] #"/staging/leuven/stg_00092/IBP_PRSproject/1_Intestinal_Bowel_Disease/"

plink_dir <- snakemake@params[[2]] #"output_data/11_test_plink/"
prsice_dir <- snakemake@params[[3]]
lasso_dir <- snakemake@params[[4]]
ldpred_dir <- snakemake@params[[5]]

target_prefix <- snakemake@params[[6]]
external_prefix <- snakemake@params[[7]]
pval_thr <- snakemake@params[[8]] #p-values used in PLINK and PRSice

range_list <- snakemake@input[[1]]
phenotype_file <- snakemake@input[[2]]

# =========================================================================================================================================================================================

setwd(project_dir)

plink_thr <- read.table(range_list, header=F)[,1]

target_plink_files <- c()
external_plink_files <- c()
for (i in 1:length(plink_thr)) {
    target_plink_files[i] <- paste(target_prefix,'.', plink_thr[i], '.profile', sep="")
    external_plink_files[i] <- paste(external_prefix,'.', plink_thr[i], '.profile', sep="")
}

# =========================================================================================================================================================================================

no_plink_thr = length(plink_thr)
no_prsice_thr = (str_count(pval_thr, ',')) + 1
# no_lasso_thr = 
# no_ldpred_thr =  

# =========================================================================================================================================================================================

setwd(project_dir)

# Phenotypes file
phenotype <- read.table(phenotype_file, header=T)
colnames(phenotype) <- c("FID", "IID","pheno") 

# PC file
pcs_file = paste('data/', target_prefix, '.eigenvec', sep="")
pcs <- read.table(pcs_file, header=F)
colnames(pcs) <- c("FID", "IID", paste0("PC",1:6)) 

covariate_file = paste('data/', target_prefix, '.cov', sep="")
covariate <- read.table(covariate_file, header=T)


# =========================================================================================================================================================================================
### PLINK ###
# =========================================================================================================================================================================================

# Get TARGET scores
setwd(paste(project_dir, plink_dir,'target_data/', sep=""))

target_prs.plink <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(target_prs.plink) <- c('FID', 'IID', 'pheno')
target_prs.plink <- read.table(target_plink_files[1], header = T)[,1:3]

for (i in 1:length(plink_thr)) {
    target_prs.plink_file <- read.table(target_plink_files[i], header = T)
    target_prs.plink_file <- target_prs.plink_file %>% select(1,2,6)
    score_col_name <- paste('SCORE', plink_thr[i], sep="")
    colnames(target_prs.plink_file) <- c('FID', 'IID', score_col_name)
    target_prs.plink <- merge(target_prs.plink, target_prs.plink_file, by = c("FID", "IID"))    
}    

#merge with the PCs
target_prs.plink <- merge(target_prs.plink, pcs, by = c("FID", "IID"))

# Get EXTERNAL scores
setwd(paste(project_dir, plink_dir,'external_data/', sep=""))

external_prs.plink <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(external_prs.plink) <- c('FID', 'IID')
external_prs.plink <- read.table(external_plink_files[1], header = T)[,1:2]

for (i in 1:length(plink_thr)) {
    external.prs.plink_file <- read.table(external_plink_files[i], header = T)
    external.prs.plink_file <- external.prs.plink_file %>% select(1,2,6)
    score_col_name_ext <- paste('SCORE', plink_thr[i], sep="")
    colnames(external.prs.plink_file) <- c('FID', 'IID', score_col_name_ext)
    external_prs.plink <- merge(external_prs.plink, external.prs.plink_file, by = c("FID", "IID"))
}    

cat("\nPLINK target and external scores uploaded into single dataframes\n")

no_plink_thr = length(plink_thr)
start_pos_ext = 3
start_pos_target = 4

ext_scores.plink <- external_prs.plink[,start_pos_ext:(start_pos_ext+no_plink_thr-1)]

ext_mean.plink <- as.data.frame(colMeans(ext_scores.plink))
ext_sd.plink <- as.data.frame(apply(ext_scores.plink, 2, sd))

colnames(ext_mean.plink) <- c("ext_mean.plink")
colnames(ext_sd.plink) <- c("ext_sd.plink")

cat("\nPLINK external scores' mean and standard deviation calculated")
cat("\nPLINK starting standardization...")

# Before standardization

std_names <- colnames(target_prs.plink)

for (sc in c(start_pos_target:(start_pos_target+no_plink_thr-1))) {
    std_names[sc] <- paste(colnames(target_prs.plink)[sc], "_std", sep="")
}

std_prs.plink <- target_prs.plink
colnames(std_prs.plink) <- std_names

for (sc in c(start_pos_target:(start_pos_target+no_plink_thr-1))) {
    col_name <- std_names[sc]
    mean <- ext_mean.plink[(sc-start_pos_target+1),1]
    stdv <- ext_sd.plink[(sc-start_pos_target+1),1]
    std_prs.plink[,sc] <- (target_prs.plink[,sc] - mean) / stdv
}

cat("\nPLINK standardization complete")

std_prs_ph.plink <- std_prs.plink
std_prs_ph.plink$PHENO <- std_prs_ph.plink$PHENO - 1

# vector with the explained variance for each threshold 
explained_var.plink <- c()
for (i in 1:length(plink_thr)) {
    
    score_var <- paste('SCORE', plink_thr[i], '_std', sep="")
    
    log_model.plink <- glm(std_prs_ph.plink$PHENO ~ pull(std_prs_ph.plink, score_var), family = binomial(link = "logit"))
    explained_var.plink[i] <- with(summary(log_model.plink), 1 - deviance/null.deviance)
}

#obtain the highest explained variance 
best_prs.plink <- which.max(explained_var.plink)

#obtain the best prs by indexing from which model had the best explained variance
best_target.plink <- std_prs_ph.plink %>% select(1, 2, 3, (start_pos_target+best_prs.plink-1))
best_target.plink <- best_target.plink[complete.cases(best_target.plink), ]
colnames(best_target.plink) <- c('FID', 'IID', 'pheno', 'std_prs')

best_thr.plink <- plink_thr[best_prs.plink]

best_thr.all <- data.frame(matrix(ncol = 2, nrow = 0))
best_thr.all <- rbind(best_thr.all, c('PLINK',best_thr.plink))
colnames(best_thr.all) <- c("tool", "best_thr")

#recode phenotype into a factor
best_target.plink$pheno <- as.factor(best_target.plink$pheno)

cat("\nPLINK scores for best threshold saved")


#need to write table for best prs
write.table(best_target.plink, snakemake@output[[1]], col.names = TRUE, row.names = FALSE, quote = FALSE)

cat("\nPLINK finished")

# =========================================================================================================================================================================================
### PRSICE ###
# =========================================================================================================================================================================================

# Get TARGET scores
setwd(paste(project_dir, prsice_dir,'target_data/', sep=""))

target_prs.prsice <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(target_prs.prsice) <- c('FID', 'IID')
target_prs.prsice_file <- read.table(paste(target_prefix,".all_score",sep=""), header = T)
target_prs.prsice <- merge(target_prs.prsice_file[,1:2], phenotype, by=c("FID", "IID"))
target_prs.prsice <- merge(target_prs.prsice, target_prs.prsice_file, by = c("FID", "IID"))
target_prs.prsice <- merge(target_prs.prsice, pcs, by = c("FID", "IID"))

# Get EXTERNAL scores
setwd(paste(project_dir, prsice_dir,'external_data/', sep=""))
external_prs.prsice <- read.table(paste(external_prefix,".all_score",sep=""), header = T)

cat("\nPRSice target and external scores uploaded into single dataframes")

no_prsice_thr = 7
start_pos_ext = 3
start_pos_target = 4

ext_scores.prsice <- external_prs.prsice[,start_pos_ext:(start_pos_ext+no_prsice_thr-1)]

ext_mean.prsice <- as.data.frame(colMeans(ext_scores.prsice))
ext_sd.prsice <- as.data.frame(apply(ext_scores.prsice, 2, sd))

colnames(ext_mean.prsice) <- c("ext_mean.prsice")
colnames(ext_sd.prsice) <- c("ext_sd.prsice")

cat("\nPRSice external scores' mean and standard deviation calculated")
cat("\nPRSice starting standardization...")

# Before standardization

std_names <- colnames(target_prs.prsice)

for (sc in c(start_pos_target:(start_pos_target+no_prsice_thr-1))) {
    std_names[sc] <- paste(colnames(target_prs.prsice)[sc], "_std", sep="")
}

std_prs.prsice <- target_prs.prsice
colnames(std_prs.prsice) <- std_names

for (sc in c(start_pos_target:(start_pos_target+no_prsice_thr-1))) {
    col_name <- std_names[sc]
    mean <- ext_mean.prsice[(sc-start_pos_target+1),1]
    stdv <- ext_sd.prsice[(sc-start_pos_target+1),1]
    std_prs.prsice[,sc] <- (target_prs.prsice[,sc] - mean) / stdv
}

cat("\nPRSice standardization complete")

std_prs_ph.prsice <- std_prs.prsice
std_prs_ph.prsice$pheno <- std_prs_ph.prsice$pheno - 1

explained_var.prsice <- c()
for (i in 1:no_prsice_thr) {    
    score_var <- colnames(std_prs_ph.prsice)[start_pos_target+i-1]
    
    log_model.prsice <- glm(std_prs_ph.prsice$pheno ~ pull(std_prs_ph.prsice, score_var), family = binomial(link = "logit"))
    explained_var.prsice[i] <- with(summary(log_model.prsice), 1 - deviance/null.deviance)
}

#obtain the highest explained variance 
best_prs.prsice <- which.max(explained_var.prsice)

#obtain the best prs by indexing from which model had the best explained variance
best_target.prsice <- std_prs_ph.prsice %>% select(1, 2, 3, (start_pos_target+best_prs.prsice-1))
best_target.prsice <- best_target.prsice[complete.cases(best_target.prsice), ]
colnames(best_target.prsice) <- c('FID', 'IID', 'pheno', 'std_prs')

# best threshold
best_thr.prsice <- colnames(target_prs.prsice)[start_pos_target+best_prs.prsice-1]
best_thr.prsice <- str_sub(best_thr.prsice, 4)
best_thr.all <- rbind(best_thr.all, c('PRSice',best_thr.prsice))
colnames(best_thr.all) <- c("tool", "best_thr")

#recode phenotype into a factor
best_target.prsice$pheno <- as.factor(best_target.prsice$pheno)

#need to write table for best prs
write.table(best_target.prsice, snakemake@output[[2]], col.names = TRUE, row.names = FALSE, quote = FALSE)

cat("\nPRSice scores for best threshold saved")
cat("\nPRSice finished")


# =========================================================================================================================================================================================
### LASSOSUM ###
# =========================================================================================================================================================================================

no_lasso_thr = 4

# Get TARGET scores
setwd(paste(project_dir,lasso_dir,'target_data/',sep=""))

target_prs.lasso <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(target_prs.lasso) <- c('FID', 'IID', 'pheno')
target_prs.lasso <- read.table(paste(target_prefix,"_prs_lasso1.txt",sep=""), header = T)[,1:3]

for (i in 1:no_lasso_thr) {
    lasso_file_name = paste(target_prefix, '_prs_lasso', i, '.txt', sep="")
    target_prs.lasso_file <- read.table(lasso_file_name, header = T)
    target_prs.lasso_file <- target_prs.lasso_file %>% select(1,2,5)
    score_col_name <- paste('best.pgs', i, sep="")
    colnames(target_prs.lasso_file) <- c('FID', 'IID', score_col_name)
    target_prs.lasso <- merge(target_prs.lasso, target_prs.lasso_file, by = c("FID", "IID"))    
}

target_prs.lasso <- merge(target_prs.lasso, pcs, by = c("FID", "IID"))

# Get EXTERNAL scores
setwd(paste(project_dir,lasso_dir,'external_data/',sep=""))

external_prs.lasso <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(external_prs.lasso) <- c('FID', 'IID')
external_prs.lasso <- read.table(paste(external_prefix,"_prs_lasso1_external.txt",sep=""), header = T)[,1:2]

for (i in 1:no_lasso_thr) {
    lasso_file_name = paste(external_prefix, '_prs_lasso', i, '_external.txt', sep="")
    external_prs.lasso_file <- read.table(lasso_file_name, header = T)
    external_prs.lasso_file <- external_prs.lasso_file %>% select(1,2,4)
    score_col_name_ext <- paste('best.pgs', i, sep="")
    colnames(external_prs.lasso_file) <- c('FID', 'IID', score_col_name_ext)
    external_prs.lasso <- merge(external_prs.lasso, external_prs.lasso_file, by = c("FID", "IID"))
}

cat("\nlassoSum target and external scores uploaded into single dataframes")

no_lasso_thr = no_lasso_thr
start_pos_ext = 3
start_pos_target = 4

ext_scores.lasso <- external_prs.lasso[,start_pos_ext:(start_pos_ext+no_lasso_thr-1)]

ext_mean.lasso <- as.data.frame(colMeans(ext_scores.lasso))
ext_sd.lasso <- as.data.frame(apply(ext_scores.lasso, 2, sd))

colnames(ext_mean.lasso) <- c("ext_mean.prsice")
colnames(ext_sd.lasso) <- c("ext_sd.prsice")

cat("\nlassoSum external scores' mean and standard deviation calculated")
cat("\nlassoSum starting standardization...")

# Before standardization
std_names <- colnames(target_prs.lasso)

for (sc in c(start_pos_target:(start_pos_target+no_lasso_thr-1))) {
    std_names[sc] <- paste(colnames(target_prs.lasso)[sc], "_std", sep="")
}

std_prs.lasso <- target_prs.lasso
colnames(std_prs.lasso) <- std_names

for (sc in c(start_pos_target:(start_pos_target+no_lasso_thr-1))) {
    col_name <- std_names[sc]
    mean <- ext_mean.lasso[(sc-start_pos_target+1),1]
    stdv <- ext_sd.lasso[(sc-start_pos_target+1),1]
    std_prs.lasso[,sc] <- (target_prs.lasso[,sc] - mean) / stdv
}

cat("\nlassoSum standardization complete")

std_prs_ph.lasso <- std_prs.lasso
std_prs_ph.lasso$pheno <- std_prs_ph.lasso$pheno - 1


explained_var.lasso <- c()
for (i in 1:no_lasso_thr) {
    
    score_var <- colnames(std_prs_ph.lasso)[start_pos_target+i-1]
    
    log_model.lasso <- glm(std_prs_ph.lasso$pheno ~ pull(std_prs_ph.lasso, score_var), family = binomial(link = "logit"))
    explained_var.lasso[i] <- with(summary(log_model.lasso), 1 - deviance/null.deviance)
}

#obtain the highest explained variance 
best_prs.lasso <- which.max(explained_var.lasso)

#obtain the best prs by indexing from which model had the best explained variance
best_target.lasso <- std_prs_ph.lasso %>% select(1, 2, 3, (start_pos_target+best_prs.lasso-1))
best_target.lasso <- best_target.lasso[complete.cases(best_target.lasso), ]
colnames(best_target.lasso) <- c('FID', 'IID', 'pheno', 'std_prs')

# best threshold
best_thr.lasso <- colnames(target_prs.lasso)[(start_pos_target+best_prs.lasso-1)]
best_thr.all <- rbind(best_thr.all, c('lassosum',best_thr.lasso))
colnames(best_thr.all) <- c("tool", "best_thr")


#recode phenotype into a factor
best_target.lasso$pheno <- as.factor(best_target.lasso$pheno)

#need to write table for best prs
write.table(best_target.lasso, snakemake@output[[3]], col.names = TRUE, row.names = FALSE, quote = FALSE)

cat("\nlassoSum scores for best threshold saved")
cat("\nlassoSum finished")



# =========================================================================================================================================================================================
### LDPRED ###
# =========================================================================================================================================================================================

## *

# Add here the code for LDpred

## *


#need to write table for best prs
write.table(best_target.ldpred, snakemake@output[[4]], col.names = TRUE, row.names = FALSE, quote = FALSE)


# =========================================================================================================================================================================================
# write table with all the thresholds
write.table(best_thr.all, snakemake@output[[5]], col.names = TRUE, row.names = FALSE, quote = FALSE)

# ===========
warnings()

