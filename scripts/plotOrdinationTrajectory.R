#################################################################################################################
# This function plots an annual trajectory over an NMDS
#################################################################################################################

 plotOrdinationTrajectory <- function(hms_prey, prey_taxonomy, startYrTS, endYrTS, preyCutoff, familyAggTS = "yes", 
                                     spp, paletteSpecies, savePlot, saveDF) {
  # Get the accum matrix
  source("./scripts/getSpMatrix.R")
  accumTS <- getSpMatrix(hms_prey = hms_prey, prey_taxonomy = prey_taxonomy, startYr = startYrTS, endYr = endYrTS,
                         endYrAc = endYrTS, preyCutoff = preyCutoff, saveFile = "no", familyAgg = familyAggTS)
  # Select specified predators
  accumTS <- subset(accumTS, Predator_Species %in% spp)
  
  # Subset so that predator/years with < 3 (or other cutoff) guts available are not included
  accumTS <- subset(accumTS, gutsPerSpPerYr >= 3)

  # Build the NMDS, this time aggregating by year
  source("./scripts/runNMDS.R")
  nmdsTS <- runNMDS(accum = accumTS, yrAgg = "yes", k = 4, try = 10, trymax = 5, saveNMDS = "no", 
                     preyCutoff = preyCutoff, familyAgg = familyAggTS, regionAgg = "yes")
  
  # Now plot, including trajectories
  data.scores <- as.data.frame(scores(nmdsTS$nmds)[[1]])
  # Add the predator name, region, year
  data.scores$grp <- nmdsTS$accumSub$Predator_Species
  data.scores$region <- nmdsTS$accumSub$Region
  data.scores$year <- nmdsTS$accumSub$Year
  # Sort
  data.scores <- data.scores[order(data.scores$grp, data.scores$year),]
  data.scores$grp <- reorder.factor(data.scores$grp, new.order=spp)
  
  # Define which predators are present here
  paletteSpecies <- subset(paletteSpecies, grp %in% data.scores$grp) 
  # Order predators by schoeners clusters
  schoeners_order <- c("Common Thresher Shark","Long-Beaked Common Dolphin",
                       "Bluefin tuna", "Albacore",
                       "Blue Shark", "Broadbill Swordfish","Bigeye Thresher Shark", "Shortfin Mako Shark",
                       "Short-Beaked Common Dolphin","Northern Right Whale Dolphin")
  paletteSpecies$grp <- reorder.factor(paletteSpecies$grp, new.order=schoeners_order)
  paletteSpecies <-paletteSpecies[order(paletteSpecies$grp),]
  # Generate an ordered color palette
  paletteSp <- paletteSpecies$col # Groups predators by color 
  data.scores$grp <- reorder.factor(data.scores$grp, new.order=schoeners_order)
  
  # Make a new DF which adds NA for years with no/insufficient samples
  missing <- data.frame(expand.grid("NMDS1" = NA, "NMDS2" = NA, "NMDS3" = NA, "NMDS4" = NA, 
                                    grp = spp, year = seq(min(accumTS$Year), max(accumTS$Year))))
  
  missing <- anti_join(missing, data.scores, by = c("grp", "year")) # Species/years with no diet info in data.scores
  data.scores.pad <- rbind(data.scores, missing)
  # Reorder
  data.scores.pad <- arrange(data.scores.pad, grp, year)
  
  # Get the last year of good data to index where the arrow head will go in plot
  maxYr <- aggregate(year ~ grp, data.scores, FUN = max, na.rm = TRUE)
  maxYr$lastYr <- maxYr$year
  data.scores <- left_join(data.scores, maxYr, by = c("grp", "year"))
  data.scores$lastYr <- ifelse(is.na(data.scores$lastYr), "no", "yes")
  
  # Generate temporal and species subsets for trajectory plots
  # Subset only including data 1998-2015
  data.scores2 <- subset(data.scores, year < 2016)
  data.scores.pad2 <- subset(data.scores.pad, year < 2016)
  
  # Subset only including swordfish and bluefin data prior up to 2015
  data.scores3 <- subset(data.scores, year < 2016 & (grp == "Bluefin tuna"| grp == "Broadbill Swordfish"))
  data.scores.pad3 <- subset(data.scores.pad2, grp == "Bluefin tuna"| grp == "Broadbill Swordfish")
  
  # Subset only including swordfish and bluefin data prior from 2015 and beyond 
  data.scores4 <- subset(data.scores,year > 2014 & (grp == "Bluefin tuna"| grp == "Broadbill Swordfish"))
  data.scores.pad4 <- subset(data.scores.pad,year > 2014 & (grp == "Bluefin tuna"| grp == "Broadbill Swordfish"))
  
  # Create subset palette for the swordgish/bluefin comparison
  paletteSpecies2 <- subset(paletteSpecies, grp %in% data.scores3$grp) 
  paletteSp2 <- paletteSpecies2$col # Groups predators by color 

  # Generate dataframe to plot vectors of families driving paritioning
  # Adding labels to family groups in ordination vectors that contribute > X dissim
  # Get the species matrix without ID cols
  spMat <- as.matrix(nmdsTS$accumSub[!names(nmdsTS$accumSub) %in% c("Predator_Species", "Predator_ID", "Year", "Region", 
                                                    "gutsPerSpPerYr")]) # exclude ID columns)
  # Fit factors to existing NMDS ordination
  NMDS.fam.fit <- envfit(nmdsTS$nmds, spMat, permutations = 999)
  # Extract scores
  NMDS.fam.scrs <- as.data.frame(scores(NMDS.fam.fit, display = "vectors")) # save species intrinsic values into dataframe
  NMDS.fam.scrs <- cbind(NMDS.fam.scrs, Prey_LPT = rownames(NMDS.fam.scrs)) # add species names to dataframe
  NMDS.fam.scrs <- cbind(NMDS.fam.scrs, pval = NMDS.fam.fit$vectors$pvals) #add pvalues to dataframe so can select sign. species
  NMDS.fam.scrs <- cbind(NMDS.fam.scrs, abrev = abbreviate(NMDS.fam.scrs$Prey_LPT, minlength = 3)) # abbreviate species names
  NMDS.fam.scrs <- subset(NMDS.fam.scrs, pval <= 0.05) # subset data to show species significant at e.g. 0.05
  NMDS.fam.scrs$hyp <- sqrt(NMDS.fam.scrs$NMDS1 ^ 2 + NMDS.fam.scrs$NMDS2 ^ 2)
  NMDS.fam.scrs2 <- subset(NMDS.fam.scrs, hyp >= 0.425) # Cutoff for important prey, imp on either axis.
  
  data.scores2$subgrp <- ifelse(data.scores2$grp == "Common Thresher Shark"|data.scores2$grp == "Broadbill Swordfish"|
                                  data.scores2$grp == "Short-Beaked Common Dolphin", 1,
                                ifelse(data.scores2$grp == "Albacore"|data.scores2$grp == "Bluefin tuna",2,
                                       ifelse(data.scores2$grp == "Long-Beaked Common Dolphin"|data.scores2$grp == "Blue Shark"|
                                                data.scores2$grp == "Shortfin Mako Shark", 3, "check")))
                                
  # Plot annual trajectories on scatter plot and facet by predator group
  pTS1 <- ggplot() + 
    stat_ellipse(data = data.scores2, aes(x = NMDS1, y = NMDS2, color = factor(grp), group = grp), size =2.5, alpha = 0.75)+
    geom_point(data = data.scores2, aes(x = NMDS1, y = NMDS2, color = factor(grp), group = grp), size = 2.5, alpha = 0.5) + 
    geom_path(data = data.scores2, aes(x = NMDS1, y = NMDS2, color = factor(grp)), lwd = 1) +
    scale_colour_manual("Predator Species", values = c(paletteSp)) +
    geom_hline(yintercept = 0, linetype = 5)+
    geom_vline(xintercept = 0, linetype = 5)+
    theme_bw(base_size = 14) +
    facet_wrap(~ subgrp, nrow=2)+
    theme(axis.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "none", legend.title = element_blank(),strip.background = element_blank(),
          strip.text.x = element_blank())
  pTS1 
  
# Plot annual trajectories on a single panel
   pTS2<- ggplot() + 
    geom_point(data = data.scores2,
               aes(x = NMDS1, y = NMDS2, color = factor(grp)), size = 4.5, alpha = 0.75) + # add the point markers
    scale_colour_manual("Predator Species", values = c(paletteSp)) +
    geom_hline(yintercept = 0, linetype = 5)+
    geom_vline(xintercept = 0, linetype = 5)+
    theme_bw(base_size = 18) +
    guides(col = guide_legend(ncol = 2, byrow = TRUE))+
    theme(axis.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "none", legend.title = element_blank(),strip.background = element_blank(),
          strip.text.x = element_blank()) +
    ## use "pad" df to only draw lines between centroids from adjacent (continuous) years
     # throws a warning about "missing rows" when you prevent paths from being drawn between points from non-sequential years
    geom_path(data = data.scores.pad2, aes(x = NMDS1, y = NMDS2, color = factor(grp)), lwd = 2,
              alpha = 0.75, # shows gaps where years missing
              arrow = arrow(ends = c("last"), angle = 20, type = "open", length = unit(0.4, "cm")), show.legend = F) +
     scale_x_continuous(limits = c(-1, 1))+
     scale_y_continuous(limits = c(-1, 1))
  pTS2
  
  # Plot annual trajectories on multiple panels to highlight change over time
  pTS3 <- ggplot() + 
    geom_point(data = data.scores3, 
               aes(x = NMDS1, y = NMDS2, color = factor(grp)), size = 3, alpha = 1) + 
    scale_colour_manual("Predator Species", values = c(paletteSp2)) +
    geom_hline(yintercept = 0, linetype = 5)+
    geom_vline(xintercept = 0, linetype = 5)+
    theme_bw(base_size = 16) +
    scale_x_continuous(limits = c(-1,0.85), breaks = c(-0.5, 0,0.5))+
    scale_y_continuous(limits = c(-.85,0.65), breaks = c(-0.5, 0, 0.5))+
    guides(col = guide_legend(ncol = 2, byrow = TRUE))+
    theme(axis.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "none", legend.title = element_blank(),strip.background = element_blank(),
          strip.text.x = element_blank()) +
    geom_path(data = data.scores.pad3, aes(x = NMDS1, y = NMDS2, color = factor(grp)), lwd = 1.5,linetype = "dashed",
              alpha=1, show.legend = F) # shows gaps where years missing
  pTS3
  
  pTS4 <- ggplot() + 
    geom_point(data = data.scores4, 
               aes(x = NMDS1, y = NMDS2, color = factor(grp)), size = 6, pch = 1, alpha = 1) + 
    scale_colour_manual("Predator Species", values = c(paletteSp2)) +
    geom_hline(yintercept = 0, linetype = 5)+
    geom_vline(xintercept = 0, linetype = 5)+
    theme_bw(base_size = 16) +
    scale_x_continuous(limits = c(-1,0.85), breaks = c(-0.5, 0,0.5))+
    scale_y_continuous(limits = c(-.85,0.65), breaks = c(-0.5, 0, 0.5))+
    guides(col = guide_legend(ncol = 2, byrow = TRUE))+
    theme(axis.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "none", legend.title = element_blank(),strip.background = element_blank(),
          strip.text.x = element_blank()) +
    geom_path(data = data.scores.pad4, aes(x = NMDS1, y = NMDS2, color = factor(grp)), lwd = 4,
              alpha = 1, # shows gaps where years missing
              show.legend = F)
  pTS4
  
  # Labeled vectors for annual NMDS
  pVectorLabs <- ggplot() + 
    geom_segment(data = NMDS.fam.scrs2, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2), arrow = arrow(length = unit(0.25, "cm")), 
                 colour = "grey10", lwd = 0.35) + # add vector arrows of significant species
    ggrepel::geom_label_repel(data = NMDS.fam.scrs2, aes(x = NMDS1, y = NMDS2, label = Prey_LPT), cex = 3, direction = "both",
                              segment.size = 0, alpha = 0.85, max.overlaps = 100, force = 1, force_pull = 50) +
    scale_colour_manual("Predator Species", values = c(paletteSp)) + 
    scale_x_continuous(limits = c(-1.5,1.5), breaks = c(-1, 0, 1))+
    scale_y_continuous(limits = c(-1.25,1.25), breaks = c(-1,-0.5, 0,  0.5,1))+
    geom_hline(yintercept = 0, linetype = 5)+
    geom_vline(xintercept = 0, linetype = 5)+
    theme_bw(base_size = 14) +
    theme(axis.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.ticks = element_blank())
  pVectorLabs
  
  # Unlabeled vectors for annual NMDS figure assembly
  pVectors <- ggplot() + 
    geom_segment(data = NMDS.fam.scrs2, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2), arrow = arrow(length = unit(0.45, "cm")), 
                 colour = "grey10", lwd = 2, alpha =0.65) + # add vector arrows of significant species
    scale_colour_manual("Predator Species", values = c(paletteSp)) + 
    scale_x_continuous(limits = c(-1.7,1.8), breaks = c(-1, 0, 1))+
    scale_y_continuous(limits = c(-1.55,1.45), breaks = c(-1, 0, 1))+
    geom_hline(yintercept = 0, linetype = 5)+
    geom_vline(xintercept = 0, linetype = 5)+
    theme_bw(base_size = 14) +
    theme(axis.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text = element_blank())
  pVectors

  # Save
  if(savePlot == "yes") {
    ggsave(paste0("./plots/ordinations/Annual_Trajectories_Faceted_byPred", startYrTS, "_", endYrTS, ".tiff"), 
           plot = pTS1, width = 4000, height = 3400, units = "px", dpi = 400, compression = "lzw")
    ggsave(paste0("./plots/ordinations/Annual_Trajectories_", startYrTS, "_", endYrTS,  ".tiff"), 
           plot = pTS2, width = 1900, height = 1930, units = "px", dpi = 400, compression = "lzw")
    ggsave(paste0("./plots/ordinations/Trajectories_BluefinSwordfish_07-14_", startYrTS, "_", endYrTS, ".tiff"), 
           plot = pTS3, width = 3000, height = 2000, units = "px", dpi = 400, compression = "lzw")
    ggsave(paste0("./plots/ordinations/Trajectories_BluefinSwordfish_post14_", startYrTS, "_", endYrTS, ".tiff"), 
           plot = pTS4, width = 3000, height = 2000, units = "px", dpi = 400, compression = "lzw")
    ggsave(paste0("./plots/ordinations/Annual_Trajectories_VectorLabels_", startYrTS, "_", endYrTS, ".tiff"), 
           plot = pVectorLabs, width = 2000, height = 2000, units = "px", dpi = 400, compression = "lzw")
    ggsave(paste0("./plots/ordinations/Annual_Trajectories_Vectors_", startYrTS, "_", endYrTS, ".tiff"), 
           plot = pVectors, width = 2000, height = 1720, units = "px", dpi = 400, compression = "lzw")
  }
  if(saveDF == "yes") {
    saveRDS(data.scores, paste0("./plots/dataFrames/plotData_plotNMDStrajectory_", startYrTS, "_", endYrTS, ".rds"))
  }
  # Return plots all together in a list
  out <- list("TrajectoriesByYear" = pTS1, "Trajectories" = pTS2, "BluefinSwordfish_0714" = pTS3, "BluefinSwordfish_post14" = pTS4, 
              "VectorLabs" = pVectorLabs, "Vectors" = pVectors, "nmdsTS" = nmdsTS) 
  return(out)
 }