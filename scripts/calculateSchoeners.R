#############################################################################################################
# Calculate Schoener's index between pairs of predators
# See Gaichas et al. 2023 http://dx.doi.org/10.1139/cjfas-2023-0093
#############################################################################################################

calculateSchoeners <- function(hms_prey, accum, savePlot, saveDF) {
  # For each prey species, calculate diff in mean propn between 2 predators, then sum those (abs) differences
  preds <- data.frame("sp1" = unique(hms_prey$Predator_Species))
  preds1 <- rep(preds$sp1, each = nrow(preds))
  preds2 <- rep(preds$sp1, times = nrow(preds))
  preds <- data.frame(cbind("sp1" = preds1, "sp2" = preds2)) # Does include comparisons against own
  preds$remove <- ifelse(preds$sp1 == preds$sp2, "remove", "ok")
  preds <- subset(preds, remove == "ok")
  preds$remove <- NULL
  # Build a df to take Schoener's coefficients. In longest format possible
  sch <- data.frame(matrix(ncol = 4, nrow = nrow(preds) * (ncol(accum) - 5)))  # 110 * 25
  colnames(sch) <- c("sp1", "sp2", "Prey_LPT", "diff")
  sch$sp1 <- rep(preds$sp1, times = (ncol(accum) - 5))
  sch$sp2 <- rep(preds$sp2, times = (ncol(accum) - 5))
  sch$Prey_LPT <- rep(colnames(accum)[6:ncol(accum)], each = nrow(preds))
  # Loop through all pairwise combinations
  for(i in 1:nrow(sch)) { 
    # Select correct prey
    toCalc <- accum[c(sch$Prey_LPT[i], "Predator_Species")] 
    # Select correct predators
    toCalc <- subset(toCalc, Predator_Species == sch$sp1[i] | Predator_Species == sch$sp2[i]) 
    agg <- aggregate(. ~ Predator_Species, toCalc, FUN = mean, na.rm = TRUE)
    sch$diff[i] <- abs(agg[1, 2] - agg[2, 2]) 
  }
  # Now calculate Schoener's index for each predator pair. e.g. alb vs pbf is 0.639. Possible values 0 --> 1 
  preds$scIndex <- NA
  for(j in 1:nrow(preds)) {
    toSum <- subset(sch, sp1 == preds$sp1[j] & sp2 == preds$sp2[j])
    preds$scIndex[j] <- 1 - 0.5 * sum(toSum$diff)
  }
  
  # Similarity matrix then cluster
  predsW <- data.frame(pivot_wider(preds, id_cols = sp1, names_from = sp2, values_from = scIndex))
  rownames(predsW) <- predsW$sp1 # Will allow correct setting of labels for clust
  # Distance matrix
  sdist <- suppressWarnings(dist(predsW[,-1]))
  clust <- hclust(sdist, method = "ward.D2")
  # Plot
  clust <- hclust(sdist, method = "ward.D2") %>% as.dendrogram %>%
    set("branches_k_color", k = 4, c("black", "black", "black", "black")) %>% set("branches_lwd",
                                                                                  c(3,3,3)) %>% set("labels_cex", c(0.75))
  # Plot in ggplot
  ggd1 <- as.ggdend(clust)
  px <- ggplot(ggd1, horiz = TRUE, offset_labels = -0.025)  + theme_bw()+ 
    theme(panel.grid = element_blank())
  px 
  if(savePlot == "yes") {
    # Save
    ggsave("./plots/SchoenersClusters.tiff", 
           plot = px, width = 2000, height = 3500, units = "px", dpi = 400, compression = "lzw")
  }
  if(saveDF == "yes") {
    saveRDS(predsW, "./plots/dataFrames/plotData_calculateSchoeners.rds")
  }
  # Return the distance matrix and the plot (which needs external editing of labels)
  out <- list("schDistMatrix" = sdist, "schCluster" = px)
  return(out) 
}