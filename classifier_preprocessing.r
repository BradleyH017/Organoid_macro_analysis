##### Bradley July 2022
##### Classifier of organoid budding/cystic classification data
##### The conda environmen
##### conda create --prefix /exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/macro_organoid_results/rf_env r-base python scikit-learn scanpy r-seurat
##### To run jupyter notbook on here:
##### pip install ipykernal


########## Looking at the correlation of covariates and classification in R
# Set up
myPaths <- .libPaths()
myPaths <- c(myPaths, "/exports/igmm/eddie/dunlop-lab/BradleyH/R/x86_64-pc-linux-gnu-library/4.0.2")
myPaths <- c(myPaths[2], myPaths[1])
.libPaths(myPaths)
setwd("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/macro_organoid_results/benchmarking_data")
library(vctrs)


# Load the data
files <- list.files()
files <- files[grep("VAC", files)]
res_list <- vector("list", length=length(files))
for(i in seq_along(res_list)){
  res_list[[i]] <- read.csv(files[i])
  res_list[[i]]$image <- gsub(".csv", "", files[i])
}

# Combine
all <- do.call(rbind, res_list)
all$image <- gsub("_all_stats_annot", "", all$image)
# write.csv(all, "all_stats_annot_all.csv")

# Remove unwanted cols
bad <- c("Mean", "StdDev", "Mode", "Median", "Min", "Max", "Skew", "Kurt", "X.Area", "Slice")
all <- all[,-which(colnames(all) %in% bad)]
# Leaves us with 28 features, inc. shape annotation and image

# Run a correlation between the featurs with one another - should we remove/merge any?
ntest <- ncol(all)-2
testnames <- colnames(all)[1:ntest]
mat <- matrix(nrow=ntest, ncol=ntest)
colnames(mat) <- testnames
rownames(mat) <- testnames
pmat <- matrix(nrow=ntest, ncol=ntest)
colnames(pmat) <- testnames
rownames(pmat) <- testnames
for(r in 1:ntest){
  for(c in 1:ntest){
    cor.res <- cor.test(all[,r], all[,c])
    mat[r,c] <- cor.res$estimate
    pmat[r,c] <- cor.res$p.value
  }
}
