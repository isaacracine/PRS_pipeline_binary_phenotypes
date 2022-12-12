# Load different libraries.
library(remotes)
# remotes::install_github("privefl/bigsnpr")
# remotes::install_github("privefl/bigstatsr")

library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

library(bigstatsr)
library(magrittr)
library(fmsb)
library(data.table)
library(runonce)

## =================================================================================
b_files <- snakemake@params[[1]]
snp_col <- snakemake@params[[2]]
noneff_allele <- snakemake@params[[3]]
eff_allele <- snakemake@params[[4]]
chrom <- snakemake@params[[5]]
pos <- snakemake@params[[6]]
beta <- snakemake@params[[7]]
pval_col <- snakemake@params[[8]]
pval <- snakemake@params[[9]]

## =================================================================================

# Read in all the different files.
sumstats <- fread(snakemake@input[[1]])
covariate <- fread(snakemake@input[[2]])
pcs <- fread(snakemake@input[[3]])
phenotype <- fread(snakemake@input[[4]])

## =================================================================================

obj.bigSNP <- snp_attach(snakemake@input[[5]]) #rds file
info <- readRDS(runonce::download_file("https://ndownloader.figshare.com/files/25503788",fname = "map_hm3_ldpred2.rds"))
bfile <- snakemake@params[[1]]
p_values <- snakemake@params[[2]]
covariate$FID <- as.character(covariate$FID)
covariate$IID <- as.character(covariate$IID)
NCORES <- nb_cores()

# Filter out hapmaps and create temp-directory
sumstats <- sumstats[, c(chrom, pos, noneff_allele, eff_allele, beta, snp_col)]
names(sumstats) <- c("chr","pos","a1","a0","beta","rsid")
tmp <- tempfile(tmpdir = "tmp-data")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
corr <- NULL
ld <- NULL
fam.order <- NULL
map <- obj.bigSNP$map[-3]
names(map) <- c("chr", "rsid", "pos", "a1", "a0")
genotype <- obj.bigSNP$genotypes
CHR <- map$chr
POS <- map$pos
info_snp <- snp_match(sumstats, map,join_by_pos = FALSE)
genotype <- obj.bigSNP$genotypes
sumstats <- sumstats[, c('chr', 'pos', 'a1', 'a0', 'beta')]
colnames(pcs) <- c("FID","IID", paste0("PC",1:6))
pcs$FID <- as.character(pcs$FID)
pcs$IID <- as.character(pcs$IID)
colnames(phenotype) <-  c('FID','IID', 'pheno_trait')
phenotype$FID <- as.character(phenotype$FID)
phenotype$IID <- as.character(phenotype$IID)
phenotype$pheno_trait <- as.character(phenotype$pheno_trait)
pheno <- merge(phenotype, covariate) %>% merge(., pcs)

POS2 <- snp_asGeneticPos(CHR, POS, dir = ".")

head(POS2)

# Calculate LD 
for (chr in 1:22) {
    ind.chr <- which(info_snp$chr == chr)
    ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
    corr0 <- snp_cor(
            genotype,
            ind.col = ind.chr2,
            ncores = NCORES,
            infos.pos = POS2[ind.chr2],
            size = 3 / 1000
        )
    if (chr == 1) {
        ld <- Matrix::colSums(corr0^2)
        corr <- as_SFBM(corr0, tmp)
    } else {
        ld <- c(ld, Matrix::colSums(corr0^2))
        corr$add_columns(corr0, nrow(corr))
    }
}




fam.order <- as.data.table(obj.bigSNP$fam)
setnames(fam.order,c("family.ID", "sample.ID"),c("FID", "IID"))

# LD score regression and modelling PRS
df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
ldsc <- snp_ldsc(   ld, 
                    length(ld), 
                    chi2 = (df_beta$beta / df_beta$beta_se)^2,
                    sample_size = df_beta$n_eff, 
                    blocks = NULL)
h2_est <- ldsc[["h2"]]

# Calculate null R2
y <- pheno[fam.order, on = c("FID", "IID")]
null.model <- paste("PC", 1:6, sep = "", collapse = "+") %>%
    paste0("Pheno~Sex", .) %>%
    as.formula %>%
    glm(., data = y, family=binomial) %>%
    summary
null.r2 <- fmsb::NagelkerkeR2(null.model)

# Prepare data for grid model
grid.param <-
    expand.grid(p = signif(seq_log(pval, 1, length.out = 1), 2),
            h2 = round(h2_est * c(0.7, 1, 1.4), 4),
            sparse = c(FALSE, TRUE))

beta_grid <-
    snp_ldpred2_grid(corr, df_beta, grid.param, ncores = NCORES)

# Obtain model PRS
genotype <- obj.bigSNP$genotypes
ind.test <- 1:nrow(genotype)
pred_grid <- big_prodMat(   genotype, 
                            beta_grid, 
                            ind.col = info_snp$`_NUM_ID_`)

# Final performance
reg.formula <- paste("PC", 1:6, sep = "", collapse = "+") %>%
    paste0("Pheno~PRS+Sex+", .) %>%
    as.formula
reg.dat_1 <- y
max.r2_1 <- 0
for(i in 1:ncol(pred_grid)){
    reg.dat_1$PRS <- pred_grid[,i]
    grid.model <- lm(reg.formula, dat=reg.dat_1) %>%
        summary  
    if(max.r2_1 < grid.model$r.squared){
        max.r2_1 <- grid.model$r.squared
    }
}
result <- data.table(grid = max.r2 - null.r2, null = null.r2)

# Creating the output. 
write.table(beta_grid, quote = FALSE, row.names = FALSE, col.names = FALSE,snakemake@output[[1]])
write.table(pred_grid, quote = FALSE, row.names = FALSE, col.names = FALSE,snakemake@output[[2]])  
write.table(grid.model, quote = FALSE, row.names = FALSE, col.names = FALSE,snakemake@output[[3]])                      
write.table(result, quote = FALSE, row.names = FALSE, col.names = FALSE,snakemake@output[[4]])                      
## End script.                      
