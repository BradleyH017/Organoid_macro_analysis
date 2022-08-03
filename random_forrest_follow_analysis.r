##### Bradley July 2022
##### Analysis of the random forrest classifications
##### Following the RF_classification.ipynb

# Set up
myPaths <- .libPaths()
myPaths <- c(myPaths, "/exports/igmm/eddie/dunlop-lab/BradleyH/R/x86_64-pc-linux-gnu-library/4.0.2")
myPaths <- c(myPaths[2], myPaths[1])
.libPaths(myPaths)
setwd("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/macro_organoid_results/results/random_forrest")
library(vctrs)

# Load the results
all <- read.csv("random_forrest_classification_results.csv")
all <- all[,-c(1,2)]
all$image <- gsub(".csv", "", all$image)

# Make the per image list
images <- levels(factor(all$image))
res_list <- vector("list", length = length(image))
for(i in seq_along(images)){
  res_list[[i]] <- all[all$image == images[i],]
}
names(res_list) <- images

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
          if(wells_want %in% well_var){
            want <- want[which(well_var == wells_want)]
            want <- do.call(rbind, want)
            temp <- data.frame(table(want$predicted_pheno))
            temp$prop <- temp$Freq/sum(temp$Freq)
            if("budding" %in% temp$Var1){
              temp <- data.frame(none=temp[temp$Var == ".",]$prop,cystic = temp[temp$Var1 == "cystic",]$prop, budding=temp[temp$Var1 == "budding",]$prop)
            } else {
              temp <- data.frame(none=temp[temp$Var == ".",]$prop,cystic = temp[temp$Var1 == "cystic",]$prop)
              temp$budding = 0
            }
          }
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

# New cols to generalise by
res_df$mouse <- unlist(strsplit(rownames(res_df), "\\."))[c(T,F,F,F,F)]
res_df$replicate <- unlist(strsplit(rownames(res_df), "\\."))[c(F,T,F,F,F)]
res_df$condition <- unlist(strsplit(rownames(res_df), "\\."))[c(F,F,T,F,F)]
res_df$day <- unlist(strsplit(rownames(res_df), "\\."))[c(F,F,F,T,F)]
res_df$well <- unlist(strsplit(rownames(res_df), "\\."))[c(F,F,F,F,T)]

# Divide into conditions - this is how we want to facet plots ultimately
conds_res <- vector("list", length = length(conds))
names(conds_res) <- conds
for(c in seq_along(conds)){
  conds_res[[c]] <- res_df[res_df$condition == conds[c],]
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
#
fold_change <- F

#
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
ppi=200
png(file="Budding_proportion_line_graphs.png", width=15*ppi, height=7.5*ppi, res=ppi)
ggplot(cond_norm_all_gen, aes(x=day, y=budding, group=mouse.rep, color=mouse)) +
    geom_point() +
    geom_line(size = 0.75) +
    facet_grid(cols = vars(condition)) +
    theme_bw() +
    ylab("Proportion of budding organoids")
dev.off()




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
write.csv(stat_all, "Budding_pvalues_within_condition.csv")


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
write.csv(stat_all, "Budding_pvalues_within_mouse.csv")

# What about the cystic
# Plot this (note: this is still un normalised so is not a fold change of this)
png(file="Cystic_proportion_line_graphs.png", width=15*ppi, height=7.5*ppi, res=ppi)
ggplot(cond_norm_all_gen, aes(x=day, y=cystic, group=mouse.rep, color=mouse)) +
    geom_point() +
    geom_line(size = 0.75) +
    facet_grid(cols = vars(condition)) +
    theme_bw() +
    ylab("Proportion of cystic organoids")
dev.off()

    library(pracma)
    cond_norm_all_gen$cond.mouse.well <- paste(cond_norm_all_gen$condition, cond_norm_all_gen$mouse.rep, sep = ".")
    mwr <- levels(factor(cond_norm_all_gen$cond.mouse.well))
    dat_list <- vector("list", length = length(mwr))
    res_df <- data.frame(cond.mouse.well = mwr, auc=rep(0,length(mwr)))
    for(r in seq_along(mwr)){
      dat_list[[r]] <- cond_norm_all_gen[cond_norm_all_gen$cond.mouse.well == mwr[r],]
      res_df[,2][r] <- trapz(dat_list[[r]]$day, dat_list[[r]]$cystic)
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
write.csv(stat_all, "Cystic_pvalues_within_condition.csv")



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
write.csv(stat_all, "Cystic_pvalues_within_mouse.csv")
