##### Bradley July 2022
##### Quick check of correlation between my own and macro budding/cystic results

# Set up
myPaths <- .libPaths()
myPaths <- c(myPaths, "/exports/igmm/eddie/dunlop-lab/BradleyH/R/x86_64-pc-linux-gnu-library/4.0.2")
myPaths <- c(myPaths[2], myPaths[1])
.libPaths(myPaths)
setwd("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/macro_organoid_results/benchmarking_data")
library(vctrs)

# Generate variable for each cystic and Budding
BH_cys <- c(62.4, 83.8, 90.6, 94.8, 59.1)
BH_bud <- c(37.6, 16.2, 9.4, 5.2, 40.9)
macro_cys <- c(63, 63.7, 73.3, 94.2, 74.4)
macro_bud <- c(37, 37.3, 26.7, 5.8, 25.6)

# calculate the correlation
cor_cys <- cor.test(BH_cys, macro_cys)
cor_bud <- cor.test(BH_bud, macro_bud)

# Calculate the lm
cys_reg <- lm(macro_cys ~ BH_cys)
bud_reg <- lm(macro_bud ~ BH_bud)

# Plot these and annotate
par(mfrow=c(1,2))
plot(BH_cys, macro_cys, main=paste0("cor=", signif(cor_cys$estimate,2), ". p=", signif(cor_cys$p.value, 2)),
  xlab="BH cystic (%)", ylab = "macro cystic (%)", pch=16, xlim=c(0,100), ylim=c(0,100))
abline(coef=c(0,1), lty="dashed", col="red")
abline(cys_reg, lty="dashed", col="blue")
plot(BH_bud, macro_bud, main=paste0("cor=", signif(cor_bud$estimate,2), ". p=", signif(cor_bud$p.value, 2)),
  xlab="BH budding (%)", ylab = "macro budding (%)", pch=16, xlim=c(0,100), ylim=c(0,100))
abline(coef=c(0,1), lty="dashed", col="red")
abline(bud_reg, lty="dashed", col="blue")

#### Finding other variables associated with the cystic/budding phenotypes
annot_files <- list.files()
annot_files <- annot_files[grep("all_stats_annot.csv", annot_files)]
res_all <- vector("list", length = length(annot_files))
for(r in seq_along(res_all)){
  res_all[[r]] <- read.csv(annot_files[r], row.names=1)
}
res <- do.call(rbind, res_all)


# Calulate the paramter-wise coefficient for association with cystic/budding
bad <- c("Mean", "StdDev", "Mode", "Min", "Max", "Median", "Skew", "Kurt", "X.Area", "Slice")
bad_cols <- which(colnames(res) %in% bad)
res_call <- res[,-bad_cols]
res_call <- res_call[res_call$Shape %in% c("cystic", "budding"),]
reg_list <- vector("list", length = ncol(res_call)[-2])
res_call$bin <- ifelse(res_call$Shape == "cystic", 0, 1)
bin_var <- res_call$bin
for(v in 1:(ncol(res_call)-2)){
  temp <- res_call[,v]
  reg_list[[v]] <- lm(temp ~ bin_var)
}

library(broom)

# Plot each
par(mfrow=c(4,7))
plist <- vector("list", length=length(reg_list))
for(v in seq_along(plist)){
  temp <- res_call[,v]
  plot(bin_var, res_call[,v], xlab="", main=colnames(res_call[v]), ylab = "")
  restemp <- tidy(reg_list[[v]])
  if(restemp$p.value[2] < 0.05){
    abline(lm(temp ~ bin_var), lty="dashed", col="red")
  }
}

# What is the correlation between significantly associated covariates
res_p <- sapply(reg_list, tidy)
assoc <- rep(F, length(res_p))
for(v in seq_along(res_p)){
  restemp <- res_p[[v]]
  if(restemp$p.value[2] < 0.05){
    assoc[v] <- T
  }
}
res_sig <- res_call[,assoc]
res_sig <- res_sig[,-which(colnames(res_sig) == "bin")]

# Matrix for results
cor_mat <- matrix(nrow=ncol(res_sig), ncol=ncol(res_sig))
dimnames(cor_mat) <- list(colnames(res_sig), colnames(res_sig))
for(r in 1:ncol(res_sig)){
  for(c in 1:ncol(res_sig)){
    res <- cor.test(res_sig[,r], res_sig[,c])
    if(res$p.value < 0.05 & abs(res$estimate > 0.5)){
      cor_mat[r,c] <- res$estimate
    } else {
      cor_mat[r,c] <- 0
    }
  }
}
write.csv(cor_mat, "pairwise_cor_assoc_vars.csv")

# Using all of the variables together
library("limma")
design_bin <- model.matrix(~bin_var)
res_test <- res_call[,-c(which(colnames(res_call) == "Shape"), which(colnames(res_call) == "bin"))]
fit_bin <- lmFit(t(res_test), design = design_bin)
fit_bin <- eBayes(fit_bin)
toptab_bin <- topTable(fit_bin, coef = 2, number = nrow(fit_bin))
toptab_bin <- toptab_bin[toptab_bin$adj.P.Val < 0.05,]
orderEffSize <- rev(order(abs(toptab_bin$logFC))) # order by effect size (absolute log-fold change)
toptab_bin[orderEffSize, ]
write.csv(toptab_bin, "toptab_bayes_lm_all_features.csv")

############### BEST
######## Trying to test a linear regression model using training and validation sets
library("minfi")
set.seed(1234)
train_prop <- floor(0.6*nrow(res_test))
train_ind <- sample(nrow(res_test), train_prop)
train_var <- res_test[train_ind,]
train_bin <- bin_var[train_ind]
valid_var <- res_test[-train_ind,]
valid_bin <- bin_var[-train_ind]

# Reducing the variables by checking for independence
cor_mat <- matrix(nrow=ncol(train_var), ncol=ncol(train_var))
dimnames(cor_mat) <- list(colnames(train_var), colnames(train_var))
for(r in 1:ncol(train_var)){
  for(c in 1:ncol(train_var)){
    res <- cor.test(train_var[,r], train_var[,c])
#    if(res$p.value < 0.05){
      cor_mat[r,c] <- res$estimate
#    } else {
#      cor_mat[r,c] <- 0
    }
  }
#}

library(ComplexHeatmap)
Heatmap(cor_mat, name = "Covariate correlations", cluster_rows=T, cluster_columns=T)

# Heatmap shows there are 6 modules of covariates
modules <- list(c("MinFeret", "Minor", "Perim.", "Width", "Height", "Area", "RawIntDen", "IntDen", "Major", "Feret"),
                c("Y", "YM", "FeretY", "BY"),
                c("X", "XM", "FeretX", "BX"),
                c("Angle", "FeretAngle"),
                c("AR"),
                c("Solidity", "Round", "Circ.")
              )
# Identify which of each module is most associated with the outcome
best_features <- rep("", length(modules))
for(m in seq_along(modules)){
  temp <- toptab_bin[rownames(toptab_bin) %in% modules[[m]],]
  temp <- temp[temp$adj.P.Val < 0.05,]
  temp <- temp[order(temp$logFC),]
  best_features[m] <- rownames(temp)[1]
}
best_features <- best_features[!is.na(best_features)]
best_features
# Shows we should only be using "Minor" "AR"    "Circ."

# Subset the training variables for these features only
train_var <- train_var[,colnames(train_var) %in% best_features]

# Generat model on training set
fit_train <- lm(train_bin ~ ., data = as.data.frame(train_var))
mean(residuals(fit_train)^2)
# 0.08311472

# Save the model
save(fit_train, file="Linear_regression_model_training.Rds")

# Check the model - Check the validation of assumptions
par(mfrow=c(2,2))
plot(fit_train)



# Now try it on the validation
mse <- function(true, prediction) {
    mean((true - prediction)^2)
}
pred_lm <- predict(fit_train, newdata = as.data.frame(valid_var))
err_lm <- mse(valid_bin, pred_lm)
err_lm
# 0.1124185

# Save the model coefficients
coefs <- coefficients(fit_train)
coefs_res <- data.frame(var=names(coefs), coefficient=coefs)
write.csv(coefs_res, "lm_pred_model_coefs.csv")

# If we round to the nearest number - how close are we for the validation set?
pred_res <- round(pred_lm)
correct <- rep(F, length(pred_res))
for(org in seq_along(pred_res)){
  if(valid_bin[org] == pred_res[org]){
    correct[org] <- T
  }
}
sum(correct)/length(correct)
# 0.8358209
# So is correct 84% of the time

# Testing on completely new data (with bad points left in)
unseen <- read.csv("VAC123_seediff3_day3.3_3_1_test_lm.csv")
unseen <- unseen[,colnames(unseen) %in% colnames(res_test)]
pred_new <- predict(fit_train, newdata = as.data.frame(unseen))
predict_round <- round(pred_new)

# And another one
unseen <- read.csv("VAC123_g2d5_day7_1_2_test_lm.csv")
unseen <- unseen[,colnames(unseen) %in% colnames(res_test)]
pred_new <- predict(fit_train, newdata = as.data.frame(unseen))
predict_round <- round(pred_new)
predict_res <- data.frame(organoid=seq(1:length(predict_round)), prediction=predict_round)
write.csv(predict_res, "VAC123_g2d5_day7_1_2_prediction_res.csv")


#~~~ Running the model on ALL Organoids. SEEM TO BE MISSING SOME SEEDIFF IMAGES
setwd("../data_all_vars")
results <- list.files()
res_list <- vector("list", length = length(results))
for(i in seq_along(results)){
  res_list[[i]] <- read.csv(results[i])
  res_list[[i]]$sample <- rep(gsub(".csv", "", results[i]), nrow(res_list[[i]]))
}
names(res_list) <- gsub(".csv", "", results)

# Remove the bad wells
count_well <- sapply(res_list, nrow)
means_well <- sapply(res_list, function(x){
  temp <- x[,which(colnames(x) == "Area")]
  mean(temp)
})
bad <- names(res_list[which(means_well > 9e4)])
remove <- which(names(res_list) %in% bad)
if(length(remove) > 0){
  res_list <- res_list[-remove]
}
highn <- names(count_well[count_well > 50])
lown <- names(count_well[count_well < 50])
if(length(lown) > 0){
  res_list <- res_list[names(res_list) %in% highn]
}

# Combine to predict
res_all <- do.call(rbind, res_list)
res_all <- res_all[,colnames(res_all) %in% colnames(res_test)]

# Use the same parameters as the size analysis to filter huge organoids. See 'organoid_size_quant_post_macro.r'
trans <- 0.4202^2
res_all$Area <- res_all$Area*(0.4202^2)
res_all <- res_all[res_all$Area < 15000,]

# Predict the Organoids
res_all$predicted_pheno <- predict(fit_train, newdata = as.data.frame(res_all))

# Remove those with absurdly high or low values. We should only be rounding to 0 and 1
res_all <- res_all[res_all$predicted_pheno < 1.5 & res_all$predicted_pheno > -0.5,]

# Plot a histogram of the actual values
hist(res_all$predicted_pheno, xlab="Phenotype prediction score", main="Distribution of prediction results")
abline(v=0.5, col="red", lty="dashed")

# Round these in a new column
res_all$predicted_bin <- round(res_all$predicted_pheno)

# Also adjust add and transform these values in the list object
for(i in seq_along(res_list)){
  res_list[[i]]$predicted_pheno <-  predict(fit_train, newdata = as.data.frame(res_list[[i]]))
  res_list[[i]] <- res_list[[i]][res_list[[i]]$predicted_pheno < 1.5 & res_list[[i]]$predicted_pheno > -0.5,]
  res_list[[i]]$predicted_bin <- round(res_list[[i]]$predicted_pheno)
}

# Setting up the list for the anlysis within biological replicates and across technical replicates
samps <- c("VAC123", "VAC131")
tog <- vector("list", length = 2)
names(tog) <- samps
reps=c("rep1", "rep2","rep3")
conds=c("growth", "g2diff", "seediff")
days=c("day1", "day3", "day5", "day7")
wells=c("well1", "well2", "well3")

# Set up
for(s in seq_along(tog)){
  # Split each sample by replicate
  tog[[s]] <- vector("list", length=length(reps))
  names(tog[[s]]) <- reps
  for(r in seq_along(reps)){
    # Split each replicate by condition
    tog[[s]][[r]] <- vector("list", length=length(conds))
    names(tog[[s]][[r]]) <- conds
    for(c in seq_along(conds)){
      # Split each condition into day
      tog[[s]][[r]][[c]] <- vector("list", length = length(days))
      names(tog[[s]][[r]][[c]]) <- days
      for(d in seq_along(days)){
        # Split each day into wells
        tog[[s]][[r]][[c]][[d]] <- vector("list", length=length(wells))
        names(tog[[s]][[r]][[c]][[d]]) <- wells
      }
    }
  }
}

# Fill list
for(s in seq_along(samps)){
  for(r in seq_along(reps)){
    for(c in seq_along(conds)){
      for(d in seq_along(days)){
        for(w in seq_along(wells)){
#          print(paste(s,r,c,d,w))
          # Mouse
          want <- res_list[grep(samps[s], names(res_list))]
          # Rep
          if(reps[r] == "rep2"){
            want <- want[grep("\\.2", names(want))]
          }
          if(reps[r] == "rep3"){
            want <- want[grep("\\.3", names(want))]
          }
          if(reps[r] == "rep1"){
            want <- want[!grepl(c("\\.2|\\.3"), names(want))]
          }
          # Conds
          if(conds[c] == "growth"){
            want <- want[grep("_g_|_g1_|_g3_|_g5_|_g7_", names(want))]
          }
          if(conds[c] == "g2diff"){
            want <- want[grep("g2d|g1d", names(want))]
          }
          if(conds[c] == "seediff"){
            want <- want[grep("seed", names(want))]
          }
          # Day
          want <- want[grep(days[d], names(want))]
          # Well
          well_var <- unlist(strsplit(names(want), "\\_"))[c(F,F,F,T,F)]
          wells_want <- gsub("well", "", wells[w])
          want <- want[which(well_var == wells_want)]
          want <- do.call(rbind, want)
          temp <- data.frame(table(want$predicted_bin))
          temp$prop <- temp$Freq/sum(temp$Freq)
          temp <- data.frame(cystic = temp[temp$Var1 == "0",]$prop, budding=temp[temp$Var1 == "1",]$prop)
          tog[[s]][[r]][[c]][[d]][[w]] <- temp
        }
        tog[[s]][[r]][[c]][[d]] <- do.call(rbind, tog[[s]][[r]][[c]][[d]])
      }
      tog[[s]][[r]][[c]] <- do.call(rbind, tog[[s]][[r]][[c]])
    }
    tog[[s]][[r]] <- do.call(rbind, tog[[s]][[r]])
  }
  tog[[s]] <- do.call(rbind, tog[[s]])
}
res_df <- do.call(rbind, tog)


# Remove those missing as are not collected
res_df <- res_df[complete.cases(res_df),]
# Rename the remaining seddiff imaged
library(data.table)
for(r in seq_along(rownames(res_df))){
  if(rownames(res_df)[r] %in% rownames(res_df[grep("day7$", rownames(res_df)),])){
    rownames(res_df)[r] <- paste0(rownames(res_df)[r], ".well1")
  }
}

#for(r in seq_along(rownames(res_df))){
#  if(rownames(res_df)[r] %like% ".day73"){
#    rep = unlist(strsplit(rownames(res_df)[r], "\\."))[c(F,T,F,F)]
#    well = "well3"
#    mouse = unlist(strsplit(rownames(res_df)[r], "\\."))[c(T,F,F,F)]
#    cond = unlist(strsplit(rownames(res_df)[r], "\\."))[c(F,F,T,F)]
#    day = unlist(strsplit(rownames(res_df)[r], "\\."))[c(F,F,F,T)]
#    day = gsub("3", "", day)
#    new = paste(mouse, rep, cond, day, well, sep = ".")
#    rownames(res_df)[r] <- new
#  }
#}

# New cols to generalise by
res_df$mouse <- unlist(strsplit(rownames(res_df), "\\."))[c(T,F,F,F,F)]
res_df$replicate <- unlist(strsplit(rownames(res_df), "\\."))[c(F,T,F,F,F)]
res_df$condition <- unlist(strsplit(rownames(res_df), "\\."))[c(F,F,T,F,F)]
res_df$day <- unlist(strsplit(rownames(res_df), "\\."))[c(F,F,F,T,F)]
res_df$well <- unlist(strsplit(rownames(res_df), "\\."))[c(F,F,F,F,T)]

# Divide into conditions - this is how we want to facet plots ultimately
conds_res <- vector("list", length = length(conds))
names(conds_res) <- conditions
for(c in seq_along(conditions)){
  conds_res[[c]] <- res_df[res_df$condition == conditions[c],]
}

##### Generalising by biological replicate - not fold changing
samps <- c("VAC123", "VAC131")
tog <- vector("list", length = 2)
names(tog) <- samps
reps=c("rep1", "rep2","rep3")
conditions=c("growth", "g2diff", "seediff")
days=c("day1", "day3", "day5", "day7")
wells=c("well1", "well2", "well3")
cond_list_gen <- conds_res

fold_change <- F

for(c in seq_along(conditions)){
  cond_list_gen[[c]] <- vector("list", length = length(samps))
  for(s in seq_along(samps)){
    cond_list_gen[[c]][[s]] <- vector("list", length=length(reps))
    for(r in seq_along(reps)){
      cond_list_gen[[c]][[s]][[r]] <- vector("list", length=length(wells))
      for(w in seq_along(wells)){
        temp <- as.data.frame(conds_res[[c]][conds_res[[c]]$mouse == samps[s] & conds_res[[c]]$replicate == reps[r] & conds_res[[c]]$well == wells[w],])
        temp <- temp[order(temp$day),]
        if(fold_change == T){
          baseline <- temp$mean[1]
          baseline_count <- temp$count[1]
          temp$mean <- (temp$mean/baseline)-1
          temp$count <- (temp$count/baseline_count)-1
        }
      cond_list_gen[[c]][[s]][[r]][[w]] <- temp
      }
      # Put wells back together
      cond_list_gen[[c]][[s]][[r]] <- do.call(rbind, cond_list_gen[[c]][[s]][[r]])
      # Now divide by day, within replciate -- CHECK THE DAYS
      temp_list <- vector("list", length = length(days))
      for(d in seq_along(days)){
        temp_list[[d]] <- cond_list_gen[[c]][[s]][[r]][cond_list_gen[[c]][[s]][[r]]$day == days[d],]
        num_cols <- unlist(lapply(temp_list[[d]], is.numeric))
        temp_list[[d]] <- sapply(temp_list[[d]][,num_cols], mean)
        # Re-add rowname
        temp_list[[d]] <- as.data.frame(t(as.data.frame(temp_list[[d]])))
        rownames(temp_list[[d]]) <- paste(c(unlist(strsplit(rownames(cond_list_gen[[c]][[s]][[r]])[1], "\\."))[c(T,T,T,F,F)], days[d]), collapse=".")
      }
      cond_list_gen[[c]][[s]][[r]] <- do.call(rbind, temp_list)
    }
    cond_list_gen[[c]][[s]] <- do.call(rbind, cond_list_gen[[c]][[s]])
  }
  cond_list_gen[[c]] <- do.call(rbind, cond_list_gen[[c]])
}
cond_norm_all_gen <- do.call(rbind, cond_list_gen)
# Re-add the other covs (lost during generalisation per biological replicate)
cond_norm_all_gen$mouse <- unlist(strsplit(rownames(cond_norm_all_gen), "\\."))[c(F,T,F,F,F)]
cond_norm_all_gen$replicate <- unlist(strsplit(rownames(cond_norm_all_gen), "\\."))[c(F,F,T,F,F)]
cond_norm_all_gen$condition <- unlist(strsplit(rownames(cond_norm_all_gen), "\\."))[c(T,F,F,F,F)]
cond_norm_all_gen$mouse.rep <- paste(cond_norm_all_gen$mouse, cond_norm_all_gen$replicate, sep = ".")
cond_norm_all_gen$condition <- factor(cond_norm_all_gen$condition, levels=c("growth", "g2diff", "seediff"))

# Add day as a numeric variable
cond_norm_all_gen$day <- as.numeric(gsub("day", "", unlist(strsplit(rownames(cond_norm_all_gen), "\\."))[c(F,F,F,F,T)]))

# Plot this (note: this is still un normalised so is not a fold change of this)
library(ggplot2)
ggplot(cond_norm_all_gen, aes(x=day, y=budding, group=mouse.rep, color=mouse)) +
    geom_point() +
    geom_line(size = 0.75) +
    facet_grid(cols = vars(condition)) +
    theme_bw() +
    ylab("Proportion of budding organoids")

# Calculate the AUC
library(pracma)
cond_norm_all_gen$cond.mouse.well <- paste(cond_norm_all_gen$condition, cond_norm_all_gen$mouse.rep, sep = ".")
mwr <- levels(factor(cond_norm_all_gen$cond.mouse.well))
dat_list <- vector("list", length = length(mwr))
res_df <- data.frame(cond.mouse.well = mwr, auc=rep(0,length(mwr)))
for(r in seq_along(mwr)){
  dat_list[[r]] <- cond_norm_all_gen[cond_norm_all_gen$cond.mouse.well == mwr[r],]
  res_df[,2][r] <- trapz(dat_list[[r]]$day, dat_list[[r]]$budding)
}

# now do t-tests within condition
library(rstatix)
res_df$condition <- unlist(strsplit(res_df$cond.mouse.well, "\\."))[c(T,F,F)]
res_df$mouse <- unlist(strsplit(res_df$cond.mouse.well, "\\."))[c(F,T,F)]
conditions=c("growth", "g2diff", "seediff")
stat.list <- vector("list", length = length(conditions))
names(stat.list)<- conditions
auc_list <- vector("list", length=length(conditions))
names(auc_list) <- conditions
for(cond in seq_along(stat.list)){
  auc_list[[cond]] <- res_df[res_df$condition == conditions[cond],]
	stat.list[[cond]] <- t_test(data=auc_list[[cond]], auc ~ mouse, paired = F)
  stat.list[[cond]] <- as.data.frame(stat.list[[cond]])
  stat.list[[cond]]$condition = conditions[cond]
}
stat_all <- do.call(rbind, stat.list)
stat_all


samps <- c("VAC123", "VAC131")
stat.list <- vector("list", length = length(samps))
names(stat.list)<- samps
auc_list <- vector("list", length=length(samps))
names(auc_list) <- samps
for(s in seq_along(stat.list)){
  auc_list[[s]] <- res_df[res_df$mouse == samps[s],]
	stat.list[[s]] <- t_test(data=auc_list[[s]], auc ~ condition, paired = F)
  stat.list[[s]] <- as.data.frame(stat.list[[s]])
  stat.list[[s]]$mouse = rep(samps[s], nrow(stat.list[[s]]))
}
stat_all <- do.call(rbind, stat.list)
