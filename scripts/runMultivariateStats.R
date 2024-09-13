#############################################################################################################
# Run multivariate analyses, then reshape Permanova and Permdisp results into triangular matrices, save them 
# Lower value of nPerm can be used for testing, can be quite slow
#############################################################################################################

runMultivariateStats <- function(accumSub, nPerm = nPerm, saveOutputs) {
  # Permanova: compare groups of objects and test the null hypothesis that the centroids and dispersion of the groups as defined 
  # by measure space are equivalent for all groups
  toPerm <- accumSub
  # Subset to 1998 - 2015
  toPerm <- subset(toPerm, Year >= 1998 & Year <= 2015)
  # Get species matrix
  spPerm <- toPerm[!names(toPerm) %in% c("Predator_Species", "Predator_ID", "Year", "Region", "gutsPerSpPerYr")]
  # Get predictors matrix
  envPerm <- toPerm[names(toPerm) %in% c("Predator_Species", "Year", "Region")]
  # Run Permanova
  perma <- adonis2(spPerm ~ Predator_Species, data = envPerm, permutations = nPerm, method = "horn", by = "margin", parallel = 4)
  # perma 
  print(paste0("Permanova complete, global p-value is ", perma$`Pr(>F)`[1]))
  
  # Pairwise tests
  permaPW <- pairwise.adonis2(spPerm ~ Predator_Species, data = envPerm, permutations = nPerm, method = "horn", by = "margin",
                              parallel = 4)
  
  # permaPW 
  print("Pairwise Permanova is complete")
  
  # Then permdisp: needs the distance matrix
  dist <- vegdist(spPerm, method = "horn")
  
  ###################################################################################################
  # Calculate mean dissim within species per year
  # Convert distance object to matrix
  distmat <- as.matrix(dist)
  # Add appropriate labels to rows/columns
  lab <- paste0(toPerm$Predator_ID, "_", toPerm$Year)
  dimnames(distmat) <- list(lab1 = lab, lab2 = lab) 
  # Melt longer, add labels
  distLong <- reshape2::melt(distmat, value.name = "distance")
  labs1 <- data.frame(str_split_fixed(distLong$lab1, "_", n = 3))
  colnames(labs1) <- c("species1", "predNo1", "year1")
  labs2 <- data.frame(str_split_fixed(distLong$lab2, "_", n = 3))
  colnames(labs2) <- c("species2", "predNo2",  "year2")
  distLong <- cbind(distLong, labs1, labs2)
  # Remove same-same comparisons
  distLong$remove <- ifelse(distLong$lab1 == distLong$lab2, "remove", "ok")
  distLong <- subset(distLong, remove == "ok")
  spp <- data.frame("sp" = unique(toPerm$Predator_Species))
  spp$mean <- spp$sd <- NA
  for (i in 1:nrow(spp)) {
    toAgg <- subset(distLong, species1 == spp$sp[i] & species2 == spp$sp[i])
    # Remove mirrors
    toAgg <- toAgg[!duplicated(apply(toAgg[, 1:2], 1, function(row) paste(sort(row), collapse = ""))), ]
    # Aggregate to show annual mean dissim within species
    distAgg <- aggregate(distance ~ year1, toAgg, FUN = mean, na.rm = TRUE)
    # Add guts/yr
    gutsPerYr <- aggregate(gutsPerSpPerYr ~ Predator_Species + Year, subset(toPerm, Predator_Species == spp$sp[i]),
                           FUN = mean)
    distAgg$Year <- as.numeric(distAgg$year1)
    distAgg <- left_join(distAgg, gutsPerYr[c("Year", "gutsPerSpPerYr")], by = "Year")
    # Just 3+ guts per year
    distAgg <- subset(distAgg, gutsPerSpPerYr > 2)
    # Mean, sd, count
    spp$mean[i] <- mean(distAgg$distance)
    spp$sd[i] <- sd(distAgg$distance)
  }
  
  
  ##########################################################################################################################
  # Complete PERMDISP
  permdisp <- betadisper(dist, group = envPerm$Predator_Species)
  permdispGlobal <- permutest(permdisp) 
  # permdispGlobal
  print(paste0("Permdisp complete, global p-value is ", permdispGlobal$tab$`Pr(>F)`[1]))
  dispPW <- data.frame(TukeyHSD(permdisp)$group) # Almost everything is significantly different
  
  # Save global results as csvs so don't need to re-run
  if(saveOutputs == "yes") {
    # Global Permanova
    write.csv(perma, paste0("./data/outputs/globalPermanovaResults_", min(toPerm$Year), "_", max(toPerm$Year), ".csv"))
    # Global Permdisp
    write.csv(permdispGlobal$tab, paste0("./data/outputs/globalPermdispResults_", min(toPerm$Year), "_", max(toPerm$Year), ".csv"))
    }
  
  ##########################################################################################################################
  # Reshape pairwise Permanova (F, p) and Permdisp results into traingular matrices, also save them. Fiddly
  # First Pairwise Permanova: reshape so is not in list form
  suppressWarnings(rm(permaPW2)) 
  for(j in 2:length(names(permaPW))) {
    out <- permaPW[[j]]
    out$Comparison <- names(permaPW)[j]
    if(exists("permaPW2")) {
      permaPW2 <- rbind(permaPW2, out)
    } else {
      permaPW2 <- out
    }
  }
  
  # Convert to triangular matrices, save both F-stats and p-values
  permaSp <- data.frame(str_split_fixed(permaPW2$Comparison, "_vs_", 2))
  colnames(permaSp) <- c("Species1", "Species2")
  permaPW3 <- cbind(permaPW2, permaSp)
  # Get triangular matrix for F
  permaPW_F <- pivot_wider(subset(permaPW3, !is.na(permaPW3$F)), id_cols = "Species1", 
                           names_from = "Species2", values_from = "F")
  # Get triangular matrix for p-values
  permaPW_P <- pivot_wider(subset(permaPW3, !is.na(permaPW3$F)), id_cols = "Species1", 
                           names_from = "Species2", values_from = "Pr(>F)")
  
  #######################################################################################################################
  # Then Permdisp
  rownames(dispPW) <- gsub("-", "_", rownames(dispPW))
  rownames(dispPW) <- gsub("g_B", "g-B", rownames(dispPW))
  rownames(dispPW) <- gsub("t_B", "t-B", rownames(dispPW))
  permDispSp <- data.frame(str_split_fixed(rownames(dispPW), "_", 2))
  colnames(permDispSp) <- c("Species1", "Species2")
  dispPW2 <- cbind(dispPW, permDispSp)
  # Get triangular matrix for p-values
  permDispPW_P <- pivot_wider(dispPW2, id_cols = "Species1", 
                           names_from = "Species2", values_from = "p.adj")
  # Save
  if(saveOutputs == "yes") {
    write.csv(permaPW_F, paste0("./data/pairwisePermanovaResults_Fstats_", 
                                min(toPerm$Year), "_", max(toPerm$Year), ".csv"))
    write.csv(permaPW_P, paste0("./data/pairwisePermanovaResults_Pvalues_", 
                                min(toPerm$Year), "_", max(toPerm$Year), ".csv"))
    write.csv(permDispPW_P, paste0("./data/pairwisePermDispResults_Pvalues_", 
                                   min(toPerm$Year), "_", max(toPerm$Year), ".csv"))
    write.csv(spp, paste0("./data/meanAnnualWithinSpDissim_", 
                          min(toPerm$Year), "_", max(toPerm$Year), ".csv"))
  }
}