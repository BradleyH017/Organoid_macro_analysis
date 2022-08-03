###### Bradley June 2022
###### Testing of organoid classifications - Using 'Budding_measure_BH.ijm' and 'Cystic_measure_BH.ijm'

##### Set up
#####
myPaths <- .libPaths()
myPaths <- c(myPaths, "/exports/igmm/eddie/dunlop-lab/BradleyH/R/x86_64-pc-linux-gnu-library/4.0.2")
myPaths <- c(myPaths[2], myPaths[1])
.libPaths(myPaths)
setwd("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/macro_organoid_results/benchmarking_data")
library(vctrs)


# Load data
res <- read.csv("Macro_shaping_benchmarking_results.csv")

# Plot the correlation for each phenotype
cor_cys <- cor.test(res$BH_cystic_proportion, res$macro_cystic_proportion)
cys_reg <-lm(macro_cystic_proportion ~ BH_cystic_proportion, data=res)
cor_bud <- cor.test(res$BH_budding_proportion, res$macro_budding_proportion)
bud_reg <-lm(macro_budding_proportion ~ BH_budding_proportion, data=res)
par(mfrow=c(1,2))
# Cystic
plot(res$BH_cystic_proportion, res$macro_cystic_proportion, xlab="BH cystic (%)", ylab="macro cystic (%)", main=paste0("cor=", signif(cor_cys$estimate,2), ". pval=", signif(cor_cys$p.value,2)))
abline(coef = c(0,1), lty="dashed", col="blue")
abline(cys_reg, lty="dashed", col="red")
# Budding
plot(res$BH_budding_proportion, res$macro_budding_proportion, xlab="BH budding (%)", ylab="macro budding (%)", main=paste0("cor=", signif(cor_bud$estimate,2), ". pval=", signif(cor_bud$p.value,2)))
abline(coef = c(0,1), lty="dashed", col="blue")
abline(bud_reg, lty="dashed", col="red")
