print("Entering the R script: global_performance_plots.R")

library(dplyr)
library(ggplot2)
library(data.table)
library(ROCR)
library(pROC)
library(pscl)
library(stringr)

# =========================================================================================================================================================================================

# # GET THEM FROM SNAKEMAKE

project_dir = snakemake@params[[1]]
target_prefix = snakemake@params[[2]]

best_target.plink_file <- snakemake@input[[1]]
best_target.prsice_file <- snakemake@input[[2]]
best_target.lasso_file <- snakemake@input[[3]]

# =========================================================================================================================================================================================

setwd(project_dir)

# PC file
pcs_file = paste('data/', target_prefix, '.eigenvec', sep="")
pcs <- read.table(pcs_file, header=F)
colnames(pcs) <- c("FID", "IID", paste0("PC",1:6)) 

# =========================================================================================================================================================================================

best_target.plink <- read.table(best_target.plink_file, header=T)
best_target.prsice <- read.table(best_target.prsice_file, header=T)
best_target.lasso <- read.table(best_target.lasso_file, header=T)

# =========================================================================================================================================================================================
### BOXPLOTS ###
# =========================================================================================================================================================================================

tmp.plink <- best_target.plink
tmp.plink$tool <- "PLINK"
tmp.prsice <- best_target.prsice
tmp.prsice$tool <- "PRSice"
tmp.lasso <- best_target.lasso
tmp.lasso$tool <- "lassosum"

tmp.all <- rbind(tmp.plink, tmp.prsice, tmp.lasso)

tmp.all$pheno <- as.factor(tmp.all$pheno)

#need to write table for all the best prs
write.table(tmp.all, snakemake@output[[1]], col.names = TRUE, row.names = FALSE, quote = FALSE)

#make boxplot
pdf(snakemake@output[[2]]) 
tmp.all %>%
    arrange(std_prs) %>%
    mutate(tool = factor(tool, levels=c('PLINK', 'PRSice', 'lassosum'))) %>%
    ggplot( aes(y = std_prs, x = pheno, fill = pheno, alpha = 0.5)) +
    labs(title = "Standardized PRS for Cases and Controls", y = "Standardized PRS", fill = 'Disease Status')+
    geom_boxplot() + guides (alpha = "none") + theme_bw() + 
    theme(plot.title = element_text(size=15), axis.title.y =  element_text(size=12), axis.text.x=element_blank(),
      axis.ticks.x=element_blank()) + scale_fill_discrete(labels = c("Controls","Cases")) +
    facet_wrap(~ tool, ncol = 3) +
    geom_signif(comparisons = list(c("0", "1")), map_signif_level=TRUE)
dev.off()



# =========================================================================================================================================================================================
### Run the ROC ###
# =========================================================================================================================================================================================

plink.roc <- roc(best_target.plink$pheno, best_target.plink$std_prs, smooth = F)
prsice.roc <- roc(best_target.prsice$pheno, best_target.prsice$std_prs, smooth = F)
lasso.roc <- roc(best_target.lasso$pheno, best_target.lasso$std_prs, smooth = F)

#do the bootstrapping 
plink.roc.ci <- ci.auc(plink.roc, method = "bootstrap", boot.n = 1000, boot.stratified = TRUE)
prsice.roc.ci <- ci.auc(prsice.roc, method = "bootstrap", boot.n = 1000, boot.stratified = TRUE)
lasso.roc.ci <- ci.auc(lasso.roc, method = "bootstrap", boot.n = 1000, boot.stratified = TRUE)

#plot the roc curve
#to add more curves just add the roc output to the list
#example: `list(lasso=lasso.roc, PLINK = plink.roc)
pdf(snakemake@output[[3]]) 
ggroc(list(PLINK=plink.roc, PRSice=prsice.roc, lassosum=lasso.roc)) + 
  labs(title = "Receiver Operator Curve Across\nConstructed Logistic Models from PRS Tools ",y = "Sensitivity (True Positive Rate)", x = "Specificity (1 - False Positive Rate)", col = "Tool") +
  theme_bw() + theme(plot.title = element_text(size=15), axis.title.y = element_text(size=12)) +
  geom_abline(intercept = 1, slope = 1)
dev.off()

# =========================================================================================================================================================================================
### Performance metrics: AUC and R2 ###
# =========================================================================================================================================================================================

### AUC ###

#write a function to allow for rounding
scaleFUN <- function(x) sprintf("%.3f", x)

#create a dataframe for the auc metric 
performance_metrics_auc <- data.frame(matrix(ncol = 4, nrow = 0))

#enter the estimates
plink_perform_auc <- c('PLINK', plink.roc.ci[1],plink.roc.ci[2],plink.roc.ci[3])
prsice_perform_auc <- c("PRSice",prsice.roc.ci[1],prsice.roc.ci[2],prsice.roc.ci[3])
lasso_perform_auc <- c("lassosum",lasso.roc.ci[1],lasso.roc.ci[2],lasso.roc.ci[3])

#combine with the empty dataframe
performance_metrics_auc <- rbind(performance_metrics_auc, plink_perform_auc, prsice_perform_auc, lasso_perform_auc)
#rename the columns of the dataframe
colnames(performance_metrics_auc) <- c("tool", "lower_95", "estimate", "upper_95")

#make the dataframe columns numeric
performance_metrics_auc$lower_95 <- as.numeric(performance_metrics_auc$lower_95)
performance_metrics_auc$estimate <- as.numeric(performance_metrics_auc$estimate)
performance_metrics_auc$upper_95 <- as.numeric(performance_metrics_auc$upper_95)

write.table(performance_metrics_auc, snakemake@output[[5]], col.names = TRUE, row.names = FALSE, quote = FALSE)

#make the auc plot
pdf(snakemake@output[[4]]) 
performance_metrics_auc %>%
    arrange(estimate) %>%
    mutate(tool = factor(tool, levels=c('PLINK', 'PRSice', 'lassosum'))) %>%
    ggplot( aes(y = estimate, x = tool, fill = tool)) + geom_bar(stat = "identity", width = 0.5, size =0.5, color = "black")+
    labs(title = bquote("Bootstrapped Estimation of AUC Value Across Tools"), y = "AUC", x = '', fill = 'Tool') +
    theme_bw() + geom_text(aes(label=scaleFUN(estimate)), vjust=0, size=3.5,y = -0.02) +
    theme(plot.title = element_text(size=15), axis.title.y = element_text(size=12), axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.05) 
dev.off()

### R2 ###

#bootstrap
n <- nrow(best_target.plink)
resamples <- 1000


best_target.plink <-  merge(best_target.plink, pcs, by =c('FID','IID'))
best_target.prsice <-  merge(best_target.prsice, pcs, by =c('FID','IID'))
best_target.lasso <-  merge(best_target.lasso, pcs, by =c('FID','IID'))

bootstrap_pseudor2_plink <- sapply(1:resamples, function(j) {
  bootstraps <- sample(c(1:n), n, TRUE)
  pR2(glm(best_target.plink$pheno[bootstraps] ~ best_target.plink$std_prs[bootstraps], family = binomial("logit")))[["McFadden"]]
})

bootstrap_pseudor2_prsice <- sapply(1:resamples, function(j) {
  bootstraps <- sample(c(1:n), n, TRUE)
  pR2(glm(best_target.prsice$pheno[bootstraps] ~ best_target.prsice$std_prs[bootstraps], family = binomial("logit")))[["McFadden"]]
})

bootstrap_pseudor2_lasso <- sapply(1:resamples, function(j) {
  bootstraps <- sample(c(1:n), n, TRUE)
  pR2(glm(best_target.lasso$pheno[bootstraps] ~ best_target.lasso$std_prs[bootstraps], family = binomial("logit")))[["McFadden"]]
})


## ==================================================================================================================
## If we you want to account for covariates, try the next code:

# bootstrap_pseudor2_plink <- sapply(1:resamples, function(j) {
#   bootstraps <- sample(c(1:n), n, TRUE)
#   pR2(glm(best_target.plink$pheno[bootstraps] ~ best_target.plink$std_prs[bootstraps] + best_target.plink$PC1[bootstraps] + best_target.plink$PC2[bootstraps] + best_target.plink$PC3[bootstraps] + best_target.plink$PC4[bootstraps] + best_target.plink$PC5[bootstraps] + best_target.plink$PC6[bootstraps], family = binomial("logit")))[["McFadden"]]
# })

# bootstrap_pseudor2_prsice <- sapply(1:resamples, function(j) {
#   bootstraps <- sample(c(1:n), n, TRUE)
#   pR2(glm(best_target.prsice$pheno[bootstraps] ~ best_target.prsice$std_prs[bootstraps] + best_target.prsice$PC1[bootstraps] + best_target.prsice$PC2[bootstraps] + best_target.prsice$PC3[bootstraps] + best_target.prsice$PC4[bootstraps] + best_target.prsice$PC5[bootstraps] + best_target.prsice$PC6[bootstraps], family = binomial("logit")))[["McFadden"]]
# })

# bootstrap_pseudor2_lasso <- sapply(1:resamples, function(j) {
#   bootstraps <- sample(c(1:n), n, TRUE)
#   pR2(glm(best_target.lasso$pheno[bootstraps] ~ best_target.lasso$std_prs[bootstraps]  + best_target.lasso$PC1[bootstraps] + best_target.lasso$PC2[bootstraps] + best_target.lasso$PC3[bootstraps] + best_target.lasso$PC4[bootstraps] + best_target.lasso$PC5[bootstraps] + best_target.lasso$PC6[bootstraps], family = binomial("logit")))[["McFadden"]]
# })
## ==================================================================================================================

#save the upper and lower bounds, as the estimate
lower_95_plink <- mean(bootstrap_pseudor2_plink) - 1.96 * sd(bootstrap_pseudor2_plink)
r_sqr_est_plink <-mean(bootstrap_pseudor2_plink)
upper_95_plink <- mean(bootstrap_pseudor2_plink) + 1.96 * sd(bootstrap_pseudor2_plink)

lower_95_prsice <- mean(bootstrap_pseudor2_prsice) - 1.96 * sd(bootstrap_pseudor2_prsice)
r_sqr_est_prsice <-mean(bootstrap_pseudor2_prsice)
upper_95_prsice <- mean(bootstrap_pseudor2_prsice) + 1.96 * sd(bootstrap_pseudor2_prsice)

lower_95_lasso <- mean(bootstrap_pseudor2_lasso) - 1.96 * sd(bootstrap_pseudor2_lasso)
r_sqr_est_lasso <-mean(bootstrap_pseudor2_lasso)
upper_95_lasso <- mean(bootstrap_pseudor2_lasso) + 1.96 * sd(bootstrap_pseudor2_lasso)


#create an empty dataframe
performance_metrics_rsqr <- data.frame(matrix(ncol = 4, nrow = 0))

#enter the estimates
plink_perform_rsqr <- c('PLINK', lower_95_plink, r_sqr_est_plink ,upper_95_plink)
prsice_perform_rsqr <- c("PRSice",lower_95_prsice, r_sqr_est_prsice ,upper_95_prsice)
lasso_perform_rsqr <- c("lassosum",lower_95_lasso, r_sqr_est_lasso ,upper_95_lasso)


#combine with the empty dataframe
performance_metrics_rsqr <- rbind(performance_metrics_rsqr, plink_perform_rsqr, prsice_perform_rsqr, lasso_perform_rsqr)
#rename the columns of the dataframe
colnames(performance_metrics_rsqr) <- c("tool", "lower_95", "estimate", "upper_95")

#make the dataframe columns numeric
performance_metrics_rsqr$lower_95 <- as.numeric(performance_metrics_rsqr$lower_95)
performance_metrics_rsqr$estimate <- as.numeric(performance_metrics_rsqr$estimate)
performance_metrics_rsqr$upper_95 <- as.numeric(performance_metrics_rsqr$upper_95)

write.table(performance_metrics_rsqr, snakemake@output[[7]], col.names = TRUE, row.names = FALSE, quote = FALSE)

#make the r-sqr plot
pdf(snakemake@output[[6]]) 
performance_metrics_rsqr %>%
    arrange(estimate) %>%
    mutate(tool = factor(tool, levels=c('PLINK', 'PRSice', 'lassosum'))) %>%
    ggplot(aes(y = estimate, x = tool, fill = tool)) + geom_bar(stat = "identity", width = 0.5, size =0.5, color = "black")+
    labs(title = bquote("Bootstrapped Estimation of"~R^2~"Value Across Tools"), y = bquote(R^2~"(Percentage of Explained Variation)"), x = '', fill = 'Tool') +
    theme_bw() + geom_text(aes(label=scaleFUN(estimate)), vjust=0, size=3.5,y = -0.005) +
    theme(plot.title = element_text(size=15), axis.title.y = element_text(size=12), axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.05) 
dev.off()

# ===========
warnings()
