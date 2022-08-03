###### Bradley June 2022
###### Analysis of organoid images, size quantified using fiji macro adapted from Laura - Sizing_Laura_BH_all.ijm
###### Make sure this is happening in a wildwest node
# Scale for each image (2.5x) is 0.4202 pixels per um

##### Set up
#####
myPaths <- .libPaths()
myPaths <- c(myPaths, "/exports/igmm/eddie/dunlop-lab/BradleyH/R/x86_64-pc-linux-gnu-library/4.0.2")
myPaths <- c(myPaths[2], myPaths[1])
.libPaths(myPaths)
setwd("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/macro_organoid_results/data")
library(vctrs)

# Read in the data as a list
# Transform as images are in pixels
trans <- 0.4202^2
files=list.files()
files <- files[grep(".csv", files)]
res_list <- vector("list", length = length(files))
names(res_list) <- gsub(".csv", "", files)
for(s in seq_along(res_list)){
  res_list[[s]] <- read.csv(files[s])
  res_list[[s]] <- res_list[[s]]$Area
  res_list[[s]] <- res_list[[s]]*trans
}


#### before summarising -- Some QC.
# 1. Plotting organoid size vs number of organoids - All replicates, all conditions, all days
count_well <- sapply(res_list, length)
means_well <- sapply(res_list, mean)

# Plot histogram of each
par(mfrow=c(2,1))
hist(count_well, xlab="Number of organoids", main="Organoid number across images")
hist(means_well, xlab="Average size of organoids", main="Mean organoid size per image")
dev.off()

# Two organoids has extremely large mean size. Pick these out
bad <- names(res_list[which(means_well > 9e4)])
remove <- which(names(res_list) %in% bad)
res_list <- res_list[-remove]
# replot
count_well <- sapply(res_list, length)
means_well <- sapply(res_list, mean)
par(mfrow=c(2,1))
hist(count_well, xlab="Number of organoids", main="Organoid number across images")
hist(means_well, xlab="Average size of organoids", main="Mean organoid size per image")
dev.off()

# Plot the correlation between these
plot(count_well, means_well, xlab="Number of organoids", ylab="Mean size of organoids")

# 2. Chuck bad quality images which contain less than 50 organoids
highn <- names(count_well[count_well > 50])
lown <- names(count_well[count_well < 50])
res_list <- res_list[names(res_list) %in% highn]

# 3. Chuck bad quality images on the basis of variance
# Images which have worked very badly (for example because of shadowing) or artefacts may have very high variance
vars <- sapply(res_list, sd)
means_well <- means_well[names(means_well) %in% names(res_list)]
par(mfrow=c(1,2))
hist(vars, xlab="Standard deviation per image")
plot(means_well, vars, xlab="Mean organoid size per image", ylab="Standard deviation per image")
dev.off()

# Find the largest organoid, check this looks okay
tops <- sapply(res_list, max)
top_tops <- max(tops)
# What sample is this
names(tops)[which(tops == max(tops))]

##### Now have finished the QC, make the master list - Can do some further QC on this
# Making the master list, sub structured
# First tier is genotype
# Then Replicate
# Then Condition
# Then Well
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



# Now set up, can subset and add to the list
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
          tog[[s]][[r]][[c]][[d]][[w]] <- unlist(want)
        }
      }
    }
  }
}


# Calculate the means
tog_means <- tog
tog_vars <- tog
tog_num <- tog
for(s in seq_along(tog)){
  for(r in seq_along(tog[[s]])){
    for(c in seq_along(tog[[s]][[r]])){
      for(d in seq_along(tog[[s]][[r]][[c]])){
        for(w in seq_along(tog[[s]][[r]][[c]][[d]])){
            temp <- mean(tog[[s]][[r]][[c]][[d]][[w]])
            tog_means[[s]][[r]][[c]][[d]][[w]] <- temp
            temp2 <- sd(tog[[s]][[r]][[c]][[d]][[w]])
            tog_vars[[s]][[r]][[c]][[d]][[w]] <- temp2
            temp3 <- length(tog[[s]][[r]][[c]][[d]][[w]])
            tog_num[[s]][[r]][[c]][[d]][[w]] <- temp3
          }
        }
      }
    }
  }

# Combine
res_df <- data.frame(unlist(tog_means))
colnames(res_df) <- "mean"
res_df$sd <- unlist(tog_vars)
res_df$count <- unlist(tog_num)
# Remove those missing as are not collected
res_df <- res_df[complete.cases(res_df),]
# Rename the remaining seddiff imaged
library(data.table)
for(r in seq_along(rownames(res_df))){
  if(rownames(res_df)[r] %like% ".day73"){
    rep = unlist(strsplit(rownames(res_df)[r], "\\."))[c(F,T,F,F)]
    well = "well3"
    mouse = unlist(strsplit(rownames(res_df)[r], "\\."))[c(T,F,F,F)]
    cond = unlist(strsplit(rownames(res_df)[r], "\\."))[c(F,F,T,F)]
    day = unlist(strsplit(rownames(res_df)[r], "\\."))[c(F,F,F,T)]
    day = gsub("3", "", day)
    new = paste(mouse, rep, cond, day, well, sep = ".")
    rownames(res_df)[r] <- new
  }
}
# New cols to generalise by
res_df$mouse <- unlist(strsplit(rownames(res_df), "\\."))[c(T,F,F,F,F)]
res_df$replicate <- unlist(strsplit(rownames(res_df), "\\."))[c(F,T,F,F,F)]
res_df$condition <- unlist(strsplit(rownames(res_df), "\\."))[c(F,F,T,F,F)]
res_df$day <- unlist(strsplit(rownames(res_df), "\\."))[c(F,F,F,T,F)]
res_df$well <- unlist(strsplit(rownames(res_df), "\\."))[c(F,F,F,F,T)]


###### Further QC - per sample/replicate
# Divide the dataframe into condition
conds <- vector("list", length = length(levels(factor(res_df$condition))))
for(c in seq_along(levels(factor(res_df$condition)))){
  conds[[c]] <- res_df[res_df$condition == levels(factor(res_df$condition))[c],]
  conds[[c]]$day <- gsub("day", "", conds[[c]]$day)
  conds[[c]]$day <- as.numeric(conds[[c]]$day)
  conds[[c]]$mouse.well <- paste(conds[[c]]$mouse, conds[[c]]$well, sep = ".")
}
conds_all <- do.call(rbind, conds)

# Plot
library(ggplot2)
conds_all$mouse.well.rep <- paste(conds_all$mouse.well, conds_all$replicate, sep = ".")
conds_all$condition <- factor(conds_all$condition, levels=c("growth", "g2diff", "seediff"))
ggplot(conds_all, aes(x=day, y=mean, group=mouse.well.rep, color=mouse)) +
    geom_point() +
    geom_line(size = 0.75) +
    facet_grid(cols = vars(condition)) +
    theme_bw() +
    ylab("Mean organoid area (um^2)")

ggplot(conds_all, aes(x=day, y=count, group=mouse.well.rep, color=mouse)) +
    geom_point() +
    geom_line(size = 0.75) +
    facet_grid(cols = vars(condition)) +
    theme_bw() +
    ylab("Number of detected organoids")


# Normalise this - so is a fold change from the day 1 value
samps <- c("VAC123", "VAC131")
tog <- vector("list", length = 2)
names(tog) <- samps
reps=c("rep1", "rep2","rep3")
conditions=c("growth", "g2diff", "seediff")
days=c("day1", "day3", "day5", "day7")
wells=c("well1", "well2", "well3")
cond_list <- conds
for(c in seq_along(conditions)){
  cond_list[[c]] <- vector("list", length = length(samps))
  for(s in seq_along(samps)){
    cond_list[[c]][[s]] <- vector("list", length=length(reps))
    for(r in seq_along(reps)){
      cond_list[[c]][[s]][[r]] <- vector("list", length=length(wells))
      for(w in seq_along(wells)){
        temp <- as.data.frame(conds[[c]][conds[[c]]$mouse == samps[s] & conds[[c]]$replicate == reps[r] & conds[[c]]$well == wells[w],])
        temp <- temp[order(temp$day),]
        if(temp$day[1] != 1){
          print("The un-anchored sample is")
          print(temp)
        }
        baseline <- temp$mean[1]
        baseline_count <- temp$count[1]
        temp$mean <- (temp$mean/baseline)-1
        temp$count <- (temp$count/baseline_count)-1
        cond_list[[c]][[s]][[r]][[w]] <- temp
      }
      cond_list[[c]][[s]][[r]] <- do.call(rbind, cond_list[[c]][[s]][[r]])
    }
    cond_list[[c]][[s]] <- do.call(rbind, cond_list[[c]][[s]])
  }
  cond_list[[c]] <- do.call(rbind, cond_list[[c]])
}
cond_norm_all <- do.call(rbind, cond_list)
cond_norm_all$mouse.well.rep <- paste(cond_norm_all$mouse.well, cond_norm_all$replicate, sep = ".")
cond_norm_all$condition <- factor(cond_norm_all$condition, levels=c("growth", "g2diff", "seediff"))

# Plot the growth fold change
ggplot(cond_norm_all, aes(x=day, y=mean, group=mouse.well.rep, color=mouse)) +
    geom_point() +
    geom_line(size = 0.75) +
    facet_grid(cols = vars(condition)) +
    theme_bw() +
    ylab("Organoid area means (Fold change)")

# Plot the number of detected organoids
ggplot(cond_norm_all, aes(x=day, y=count, group=mouse.well.rep, color=mouse)) +
    geom_point() +
    geom_line(size = 0.75) +
    facet_grid(cols = vars(condition)) +
    theme_bw() +
    ylab("Organoids detected (Fold change)")


#### Generalising the results per biological replicate - So now wat to merge the
samps <- c("VAC123", "VAC131")
tog <- vector("list", length = 2)
names(tog) <- samps
reps=c("rep1", "rep2","rep3")
conditions=c("growth", "g2diff", "seediff")
days=c("day1", "day3", "day5", "day7")
wells=c("well1", "well2", "well3")
cond_list_gen <- conds
for(c in seq_along(conditions)){
  cond_list_gen[[c]] <- vector("list", length = length(samps))
  for(s in seq_along(samps)){
    cond_list_gen[[c]][[s]] <- vector("list", length=length(reps))
    for(r in seq_along(reps)){
      cond_list_gen[[c]][[s]][[r]] <- vector("list", length=length(wells))
      for(w in seq_along(wells)){
        temp <- as.data.frame(conds[[c]][conds[[c]]$mouse == samps[s] & conds[[c]]$replicate == reps[r] & conds[[c]]$well == wells[w],])
        temp <- temp[order(temp$day),]
        baseline <- temp$mean[1]
        baseline_count <- temp$count[1]
        temp$mean <- (temp$mean/baseline)-1
        temp$count <- (temp$count/baseline_count)-1
        cond_list_gen[[c]][[s]][[r]][[w]] <- temp
      }
      # Put wells back together
      cond_list_gen[[c]][[s]][[r]] <- do.call(rbind, cond_list_gen[[c]][[s]][[r]])
      # Now divide by day, within replciate
      temp_list <- vector("list", length = length(days))
      for(d in seq_along(days)){
        dnum <- gsub("day", "", days[d])
        temp_list[[d]] <- cond_list_gen[[c]][[s]][[r]][cond_list_gen[[c]][[s]][[r]]$day == dnum,]
        num_cols <- unlist(lapply(temp_list[[d]], is.numeric))
        temp_list[[d]] <- sapply(temp_list[[d]][,num_cols], mean)
        # Re-add rowname
        temp_list[[d]] <- as.data.frame(t(as.data.frame(temp_list[[d]])))
        rownames(temp_list[[d]]) <- paste(c(unlist(strsplit(rownames(cond_list_gen[[c]][[s]][[r]])[1], "\\."))[c(T,T,T,F,F)], paste0("day",temp_list[[d]]$day)), collapse=".")
      }
      cond_list_gen[[c]][[s]][[r]] <- do.call(rbind, temp_list)
    }
    cond_list_gen[[c]][[s]] <- do.call(rbind, cond_list_gen[[c]][[s]])
  }
  cond_list_gen[[c]] <- do.call(rbind, cond_list_gen[[c]])
}
cond_norm_all_gen <- do.call(rbind, cond_list_gen)
# Re-add the other covs (lost during generalisation per biological replicate)
cond_norm_all_gen$mouse <- unlist(strsplit(rownames(cond_norm_all_gen), "\\."))[c(T,F,F,F)]
cond_norm_all_gen$replicate <- unlist(strsplit(rownames(cond_norm_all_gen), "\\."))[c(F,T,F,F)]
cond_norm_all_gen$condition <- unlist(strsplit(rownames(cond_norm_all_gen), "\\."))[c(F,F,T,F)]
cond_norm_all_gen$mouse.rep <- paste(cond_norm_all_gen$mouse, cond_norm_all_gen$replicate, sep = ".")
cond_norm_all_gen$condition <- factor(cond_norm_all_gen$condition, levels=c("growth", "g2diff", "seediff"))

# Plot the growth fold change per biological replicate
ppi=300
png(file="../results/FC_growth_curves_bio_reps.png", width=12*ppi, height=6*ppi, res=ppi)
ggplot(cond_norm_all_gen, aes(x=day, y=mean, group=mouse.rep, color=mouse)) +
    geom_point() +
    geom_line(size = 0.75) +
    facet_grid(cols = vars(condition)) +
    theme_bw() +
    ylab("Organoid area means (Fold change)")
dev.off()


# Calculating the area under the curves per mouse.rep, within conditions
library(pracma)
cond_norm_all_gen$cond.mouse.well <- paste(cond_norm_all_gen$condition, cond_norm_all_gen$mouse.rep, sep = ".")
mwr <- levels(factor(cond_norm_all_gen$cond.mouse.well))
dat_list <- vector("list", length = length(mwr))
res_df <- data.frame(cond.mouse.well = mwr, auc=rep(0,length(mwr)))
for(r in seq_along(mwr)){
  dat_list[[r]] <- cond_norm_all_gen[cond_norm_all_gen$cond.mouse.well == mwr[r],]
  res_df[,2][r] <- trapz(dat_list[[r]]$day, dat_list[[r]]$mean)
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
write.csv(stat_all, "../results/t_test_within_condition_across_mouse.csv")

# Now do t-tests across condition, within mouse
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
write.csv(stat_all, "../results/t_test_across_condition_within_mouse.csv")


# Plotting all the data - individual organoids
# Need to melt it all together
res_list_temp <- res_list
for(i in seq_along(res_list_temp)){
  res_list_temp[[i]] <- data.frame(sample=rep(names(res_list)[i], length(res_list[[i]])), value=res_list[[i]])
}
melt_temp <- do.call(rbind, res_list_temp)
melt_temp$mouse <- unlist(strsplit(melt_temp$sample, "\\_"))[c(T,F,F,F,F)]
melt_temp$condition <- rep("", nrow(melt_temp))
melt_temp$day = rep("", nrow(melt_temp))
melt_temp$rep = rep("", nrow(melt_temp))
for(r in 1:nrow(melt_temp)){
  apart <- unlist(strsplit(melt_temp$sample[r], "\\_"))
  day_rep <- apart[grep("day", apart)]
  if(day_rep %in% paste("day", c("1.3", "3.3", "5.3", "7.3"), sep = "")){
    melt_temp$day[r] <- day_rep
    melt_temp$rep[r] <- "rep3"
  } else {
  if(day_rep %in% paste("day", c("1.2", "3.2", "5.2", "7.2"), sep = "")){
      melt_temp$day[r] <- day_rep
      melt_temp$rep[r] <- "rep2"
    }
  if(melt_temp$rep[r] == ""){
    melt_temp$day[r] <- day_rep
    melt_temp$rep[r] <- "rep1"
  }
  temp <- apart[grep("g|g1|g3|g5|g7|g2d|g1d|seed", apart)]
  if(temp %like% "g1d|g2d"){
    melt_temp$condition[r] <- "g2diff"
  }
  if(temp %like% "seed"){
    melt_temp$condition[r] <- "seediff"
  }
  if(melt_temp$condition[r] == ""){
    melt_temp$condition <- "growth"
    }
  }
}
melt_temp$day <- gsub("\\.2", "", melt_temp$day)
melt_temp$day <- gsub("\\.3", "", melt_temp$day)
melt_temp$mouse.rep <- paste(melt_temp$mouse, melt_temp$rep, sep = ".")
write.csv(melt_temp, "../temp/melt_organoid_vals_temp.csv")

# Violin plot
melt_temp$day <- factor(melt_temp$day, levels = c("day1", "day3", "day5", "day7"))
melt_temp$log_value <- log10(melt_temp$value)
library(ggbeeswarm)
ggplot(melt_temp, aes(x=day, y=log_value, fill=mouse)) +
  geom_violin() +
  geom_violin(alpha = 0.5) +
  geom_quasirandom(inherit.aes = TRUE) +
#  geom_boxplot(width=1, color="grey") +
#  geom_dotplot(binaxis="y", stackdir='center', dotsize=0.005) +
#  geom_jitter(shape=16, position=position_jitter(0.05), size=0.05) +
  facet_grid(cols = vars(condition)) +
  theme_bw() +
  ylab("Organoid area (um^2)")


# Beeswarm plot
melt_temp$mouse.condition <- paste(melt_temp$mouse, melt_temp$condition, sep = ".")
beecond <- levels(factor(melt_temp$mouse.condition))
beedata <- vector("list", length = length(beecond))
plots <- vector("list", length=length(beecond))
melt_temp$mouse.condition <- factor(melt_temp$mouse.condition, levels = c("VAC123.growth", "VAC123.g2diff", "VAC123.seediff", "VAC131.growth", "VAC131.g2diff", "VAC131.seediff"))
par(mfrow=c(2,3))
for(bc in seq_along(beecond)){
  beedata[[bc]] <- melt_temp[melt_temp$mouse.condition == beecond[bc],]
  beeswarm(value ~ day,
    data=beedata[[bc]],
    pch = 1, cex=0.05, do.plot=T, main=beecond[bc], ylim=c(0,max(melt_temp$value)))
}


##### Post beeswarm Analysis
# This shows we have some measurements that are exceptionally large (and likley include > 1 organoid)
# If we draw a cut off of 20,000um^2, how many organoids do we lose? - 90, largest is still incorrect
# Using 15,000
melt_temp <- read.csv("../temp/melt_organoid_vals_temp.csv", row.names=1)
old <- nrow(melt_temp)
melt_filt <- melt_temp[melt_temp$value < 15000,]
new <- nrow(melt_filt)
old-new
rownames(melt_filt[melt_filt$value == max(melt_filt$value),])

# Repeating with these Organoids
melt_filt$mouse.condition <- paste(melt_temp$mouse, melt_temp$condition, sep = ".")
beecond <- levels(factor(melt_filt$mouse.condition))
beedata <- vector("list", length = length(beecond))
plots <- vector("list", length=length(beecond))
melt_filt$mouse.condition <- factor(melt_filt$mouse.condition, levels = c("VAC123.growth", "VAC123.g2diff", "VAC123.seediff", "VAC131.growth", "VAC131.g2diff", "VAC131.seediff"))
par(mfrow=c(2,3))
for(bc in seq_along(beecond)){
  beedata[[bc]] <- melt_filt[melt_filt$mouse.condition == beecond[bc],]
  beeswarm(value ~ day,
    data=beedata[[bc]],
    pch = 1, cex=0.05, do.plot=T, main=beecond[bc], ylim=c(0,max(melt_filt$value)))
}

# try with test value - 5000
melt_test <- melt_temp[melt_temp$value < 5000,]
old <- nrow(melt_temp)
new <- nrow(melt_test)
old-new
melt_test$mouse.condition <- paste(melt_test$mouse, melt_test$condition, sep = ".")
beecond <- levels(factor(melt_test$mouse.condition))
beedata <- vector("list", length = length(beecond))
plots <- vector("list", length=length(beecond))
melt_test$mouse.condition <- factor(melt_test$mouse.condition, levels = c("VAC123.growth", "VAC123.g2diff", "VAC123.seediff", "VAC131.growth", "VAC131.g2diff", "VAC131.seediff"))
par(mfrow=c(2,3))
for(bc in seq_along(beecond)){
  beedata[[bc]] <- melt_test[melt_test$mouse.condition == beecond[bc],]
  beeswarm(value ~ day,
    data=beedata[[bc]],
    pch = 1, cex=0.05, do.plot=T, main=beecond[bc], ylim=c(0,max(melt_test$value)))
}
######################################## ~~~~~~~~~~~~~~~~~~~~~~
##### Subsetting the res_list for only the organoids < 15,000 in Area
for(i in seq_along(res_list)){
  res_list[[i]] <- res_list[[i]][res_list[[i]] < 15000]
}


# Now produce the plot
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



# Now set up, can subset and add to the list
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
          tog[[s]][[r]][[c]][[d]][[w]] <- unlist(want)
        }
      }
    }
  }
}


# Calculate the means
tog_means <- tog
tog_vars <- tog
tog_num <- tog
for(s in seq_along(tog)){
  for(r in seq_along(tog[[s]])){
    for(c in seq_along(tog[[s]][[r]])){
      for(d in seq_along(tog[[s]][[r]][[c]])){
        for(w in seq_along(tog[[s]][[r]][[c]][[d]])){
            temp <- mean(tog[[s]][[r]][[c]][[d]][[w]])
            tog_means[[s]][[r]][[c]][[d]][[w]] <- temp
            temp2 <- sd(tog[[s]][[r]][[c]][[d]][[w]])
            tog_vars[[s]][[r]][[c]][[d]][[w]] <- temp2
            temp3 <- length(tog[[s]][[r]][[c]][[d]][[w]])
            tog_num[[s]][[r]][[c]][[d]][[w]] <- temp3
          }
        }
      }
    }
  }

# Combine
res_df <- data.frame(unlist(tog_means))
colnames(res_df) <- "mean"
res_df$sd <- unlist(tog_vars)
res_df$count <- unlist(tog_num)
# Remove those missing as are not collected
res_df <- res_df[complete.cases(res_df),]
# Rename the remaining seddiff imaged
library(data.table)
for(r in seq_along(rownames(res_df))){
  if(rownames(res_df)[r] %like% ".day73"){
    rep = unlist(strsplit(rownames(res_df)[r], "\\."))[c(F,T,F,F)]
    well = "well3"
    mouse = unlist(strsplit(rownames(res_df)[r], "\\."))[c(T,F,F,F)]
    cond = unlist(strsplit(rownames(res_df)[r], "\\."))[c(F,F,T,F)]
    day = unlist(strsplit(rownames(res_df)[r], "\\."))[c(F,F,F,T)]
    day = gsub("3", "", day)
    new = paste(mouse, rep, cond, day, well, sep = ".")
    rownames(res_df)[r] <- new
  }
}
# New cols to generalise by
res_df$mouse <- unlist(strsplit(rownames(res_df), "\\."))[c(T,F,F,F,F)]
res_df$replicate <- unlist(strsplit(rownames(res_df), "\\."))[c(F,T,F,F,F)]
res_df$condition <- unlist(strsplit(rownames(res_df), "\\."))[c(F,F,T,F,F)]
res_df$day <- unlist(strsplit(rownames(res_df), "\\."))[c(F,F,F,T,F)]
res_df$well <- unlist(strsplit(rownames(res_df), "\\."))[c(F,F,F,F,T)]

conds <- vector("list", length = length(levels(factor(res_df$condition))))
for(c in seq_along(levels(factor(res_df$condition)))){
  conds[[c]] <- res_df[res_df$condition == levels(factor(res_df$condition))[c],]
  conds[[c]]$day <- gsub("day", "", conds[[c]]$day)
  conds[[c]]$day <- as.numeric(conds[[c]]$day)
  conds[[c]]$mouse.well <- paste(conds[[c]]$mouse, conds[[c]]$well, sep = ".")
}
conds_all <- do.call(rbind, conds)
samps <- c("VAC123", "VAC131")
tog <- vector("list", length = 2)
names(tog) <- samps
reps=c("rep1", "rep2","rep3")
conditions=c("growth", "g2diff", "seediff")
days=c("day1", "day3", "day5", "day7")
wells=c("well1", "well2", "well3")
cond_list_gen <- conds
for(c in seq_along(conditions)){
  cond_list_gen[[c]] <- vector("list", length = length(samps))
  for(s in seq_along(samps)){
    cond_list_gen[[c]][[s]] <- vector("list", length=length(reps))
    for(r in seq_along(reps)){
      cond_list_gen[[c]][[s]][[r]] <- vector("list", length=length(wells))
      for(w in seq_along(wells)){
        temp <- as.data.frame(conds[[c]][conds[[c]]$mouse == samps[s] & conds[[c]]$replicate == reps[r] & conds[[c]]$well == wells[w],])
        temp <- temp[order(temp$day),]
        baseline <- temp$mean[1]
        baseline_count <- temp$count[1]
        temp$mean <- (temp$mean/baseline)-1
        temp$count <- (temp$count/baseline_count)-1
        cond_list_gen[[c]][[s]][[r]][[w]] <- temp
      }
      # Put wells back together
      cond_list_gen[[c]][[s]][[r]] <- do.call(rbind, cond_list_gen[[c]][[s]][[r]])
      # Now divide by day, within replciate
      temp_list <- vector("list", length = length(days))
      for(d in seq_along(days)){
        dnum <- gsub("day", "", days[d])
        temp_list[[d]] <- cond_list_gen[[c]][[s]][[r]][cond_list_gen[[c]][[s]][[r]]$day == dnum,]
        num_cols <- unlist(lapply(temp_list[[d]], is.numeric))
        temp_list[[d]] <- sapply(temp_list[[d]][,num_cols], mean)
        # Re-add rowname
        temp_list[[d]] <- as.data.frame(t(as.data.frame(temp_list[[d]])))
        rownames(temp_list[[d]]) <- paste(c(unlist(strsplit(rownames(cond_list_gen[[c]][[s]][[r]])[1], "\\."))[c(T,T,T,F,F)], paste0("day",temp_list[[d]]$day)), collapse=".")
      }
      cond_list_gen[[c]][[s]][[r]] <- do.call(rbind, temp_list)
    }
    cond_list_gen[[c]][[s]] <- do.call(rbind, cond_list_gen[[c]][[s]])
  }
  cond_list_gen[[c]] <- do.call(rbind, cond_list_gen[[c]])
}
cond_norm_all_gen <- do.call(rbind, cond_list_gen)
# Re-add the other covs (lost during generalisation per biological replicate)
cond_norm_all_gen$mouse <- unlist(strsplit(rownames(cond_norm_all_gen), "\\."))[c(T,F,F,F)]
cond_norm_all_gen$replicate <- unlist(strsplit(rownames(cond_norm_all_gen), "\\."))[c(F,T,F,F)]
cond_norm_all_gen$condition <- unlist(strsplit(rownames(cond_norm_all_gen), "\\."))[c(F,F,T,F)]
cond_norm_all_gen$mouse.rep <- paste(cond_norm_all_gen$mouse, cond_norm_all_gen$replicate, sep = ".")
cond_norm_all_gen$condition <- factor(cond_norm_all_gen$condition, levels=c("growth", "g2diff", "seediff"))

# Plot the growth fold change per biological replicate
ppi=300
png(file="../results/FC_growth_curves_bio_reps_filt_15k.png", width=12*ppi, height=6*ppi, res=ppi)
ggplot(cond_norm_all_gen, aes(x=day, y=mean, group=mouse.rep, color=mouse)) +
    geom_point() +
    geom_line(size = 0.75) +
    facet_grid(cols = vars(condition)) +
    theme_bw() +
    ylab("Organoid area means (Fold change)")
dev.off()

# Calculate auc
library(pracma)
cond_norm_all_gen$cond.mouse.well <- paste(cond_norm_all_gen$condition, cond_norm_all_gen$mouse.rep, sep = ".")
mwr <- levels(factor(cond_norm_all_gen$cond.mouse.well))
dat_list <- vector("list", length = length(mwr))
res_df <- data.frame(cond.mouse.well = mwr, auc=rep(0,length(mwr)))
for(r in seq_along(mwr)){
  dat_list[[r]] <- cond_norm_all_gen[cond_norm_all_gen$cond.mouse.well == mwr[r],]
  res_df[,2][r] <- trapz(dat_list[[r]]$day, dat_list[[r]]$mean)
}

# Do the tests
library(rstatix)
conditions=c("growth", "g2diff", "seediff")
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
write.csv(stat_all, "../results/t_test_within_condition_across_mouse_filt15k.csv")

# Now do t-tests across condition, within mouse
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
write.csv(stat_all, "../results/t_test_across_condition_within_mouse_filt15k.csv")
