#################################################################################################################
# This function plots a scatterplot of an NMDS
# Function requires an NMDS object, and a species matrix.
#################################################################################################################

plotScatter <- function(myNMDS, myAccum, preyCutoff, familyAgg, yrAgg, paletteSpecies, 
                        savePlot, saveDF) {
  # Using the scores function from vegan to extract the site scores and convert to a data.frame
  data.scores <- as.data.frame(scores(myNMDS) [[1]])  ### I added list index, required
  # Create a column of site names, from the rownames of data.scores
  data.scores$site <- rownames(data.scores) 
  # Add the grp variable created earlier
  data.scores$grp <- myAccum$Predator_Species
  # Add region and year too
  data.scores$region <- myAccum$Region
  data.scores$year <- myAccum$Year
  
  # Optional: subset so that predator/years with < 3 (or other cutoff) guts available are not included
  data.scores$gutsPerSpPerYr <- myAccum$gutsPerSpPerYr
  
  # Define which predators are present here
  paletteSpecies <- subset(paletteSpecies, grp %in% data.scores$grp) 
  # Order predators by schoeners clusters
  schoeners_order <- c("Broadbill Swordfish", "Blue Shark","Common Thresher Shark","Long-Beaked Common Dolphin",
                       "Albacore", "Shortfin Mako Shark","Bigeye Thresher Shark","Short-Beaked Common Dolphin",
                       "Bluefin tuna",""," ", "Northern Right Whale Dolphin")
  paletteSpecies$grp <- reorder.factor(paletteSpecies$grp, new.order=schoeners_order)
  paletteSpecies <-paletteSpecies[order(paletteSpecies$grp),]
  paletteSp <- paletteSpecies$col # Groups predators by color 
  data.scores$grp <- reorder.factor(data.scores$grp, new.order=schoeners_order)
  
  # Plot a simple scatter 
  p1 <- ggplot() + 
    geom_point(data = data.scores, aes(x = NMDS1, y = NMDS2, color = factor(grp)), size = 1.5, alpha = 0.5) + # add the point markers
    stat_ellipse(data = data.scores, aes(x = NMDS1, y = NMDS2, color = factor(grp)), size =2.5, alpha = 0.75)+
    scale_colour_manual("Predator Species", values = c(paletteSp)) +
    geom_hline(yintercept = 0, linetype = 5)+
    geom_vline(xintercept = 0, linetype = 5)+
    theme_bw(base_size = 14) +
    guides(col = guide_legend(ncol = 3, byrow = TRUE))+
    theme(axis.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "bottom", legend.title = element_blank(),strip.background = element_blank(),
          strip.text.x = element_blank())
  
  # plot without legend for easy assembly in Photoshop (faceted)
  p2 <- ggplot() + 
    geom_point(data = data.scores, aes(x = NMDS1, y = NMDS2, color = factor(grp)), size = 2, alpha = 0.5) + # add the point markers
    stat_ellipse(data = data.scores, aes(x = NMDS1, y = NMDS2, color = factor(grp)), size =1.5, alpha = 1)+
    scale_colour_manual("Predator Species", values = c(paletteSp)) +
    scale_x_reverse()+
    geom_hline(yintercept = 0, linetype = 5)+
    geom_vline(xintercept = 0, linetype = 5)+
    theme_bw(base_size = 14) +
    facet_wrap(~grp, nrow = 3, ncol = 4, drop = F)+
    theme(axis.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "none", strip.background = element_blank(),
          strip.text.x = element_blank())
  p2
  
  # Get legend in case needed for figure assembly (requires ggpubr package)
  leg <- get_legend(p1)
  # Convert to a ggplot and print
  p3 <- as_ggplot(leg)

  # Get centroid and sd within species/year: try to show how much of variability is covered in each year
  data.scores.agg <- aggregate(cbind(NMDS1, NMDS2) ~ grp+year+gutsPerSpPerYr, data.scores, FUN = median, na.rm = TRUE)
  data.scores.agg2 <- subset(data.scores.agg, gutsPerSpPerYr > 2)
  data.scores.agg3 <- aggregate(cbind(NMDS1, NMDS2) ~ grp, data.scores.agg2, FUN = median, na.rm = TRUE)

  p4 <- ggplot() + 
    stat_ellipse(data = data.scores, aes(x = NMDS1, y = NMDS2, color = factor(grp), group = year),
                 size = 1, alpha = 0.75) + # Many sp/yrs will have too few points to calculate an ellipse
    geom_hline(yintercept = 0, linetype = 5)+
    geom_vline(xintercept = 0, linetype = 5)+
    geom_point(data = data.scores.agg3, aes(x = NMDS1, y = NMDS2, fill = factor(grp)), size = 4, pch = 21, 
               alpha =0.4) + # Just the overall spp mean
    scale_colour_manual("Predator Species", values = c(paletteSp)) +
    scale_fill_manual("Predator Species", values = c(paletteSp)) +
    theme_bw(base_size = 14) +
    facet_wrap(~grp, nrow = 3, ncol = 4)+
    theme(axis.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "none", strip.background = element_blank(),
          strip.text.x = element_blank())
  p4
  
  p5 <- ggplot() + 
    stat_ellipse(data = data.scores.agg2, aes(x = NMDS1, y = NMDS2, color = factor(grp)), size =1.5, alpha = 0.75)+
    geom_hline(yintercept = 0, linetype = 5)+
    geom_vline(xintercept = 0, linetype = 5)+
    geom_path(data = data.scores.agg2, aes(x = NMDS1, y = NMDS2, color = factor(grp), group =year*grp), 
              lwd = 1, # no gaps where years missing
    arrow = arrow(ends = c("last"), angle = 20, type = "open", length = unit(0.35, "cm")), show.legend = F, alpha =0.85) +
    scale_colour_manual("Predator Species", values = c(paletteSp)) +
    scale_fill_manual("Predator Species", values = c(paletteSp)) +
    theme_bw(base_size = 14) +
    facet_wrap(~grp, nrow = 3, ncol = 4)+
    theme(axis.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "none", strip.background = element_blank(),
          strip.text.x = element_blank())
  p5

    # Save
  if(savePlot == "yes") {
     ggsave(paste0("./plots/ordinations/simpleScatter_", preyCutoff, "_", "familyAgg", familyAgg, 
                   "_yrAgg", yrAgg, "_", min(myAccum$Year), "_", max(myAccum$Year),  ".tiff"), 
                   plot = p1, width = 4000, height = 2000, units = "px", dpi = 400, compression = "lzw")
    ggsave(paste0("./plots/ordinations/simpleScatter_byPredator", preyCutoff, "_", "familyAgg", familyAgg, 
                  "_yrAgg", yrAgg, "_", min(myAccum$Year), "_", max(myAccum$Year),  ".tiff"), 
           plot = p2, width = 3341, height = 2214, units = "px", dpi = 400, compression = "lzw")
    ggsave(paste0("./plots/ordinations/simpleScatter_LEGEND", preyCutoff, "_", "familyAgg", familyAgg, 
                  "_yrAgg", yrAgg, "_", min(myAccum$Year), "_", max(myAccum$Year),"_LEGEND" , ".tiff"), 
           plot = p3, width = 2400, height = 1000, units = "px", dpi = 400, compression = "lzw")
    ggsave(paste0("./plots/ordinations/scatterYearEllipses_", preyCutoff, "_", "familyAgg", familyAgg, 
                  "_yrAgg", yrAgg, "_", min(myAccum$Year), "_", max(myAccum$Year), ".tiff"), 
           plot = p4, width = 4000, height = 2000, units = "px", dpi = 400, compression = "lzw")
    ggsave(paste0("./plots/ordinations/scatterYearTrajectories_", preyCutoff, "_", "familyAgg", familyAgg, 
                  "_yrAgg", yrAgg, "_", min(myAccum$Year), "_", max(myAccum$Year),  ".tiff"), 
           plot = p5, width = 4000, height = 2000, units = "px", dpi = 400, compression = "lzw")
  }
  if(saveDF == "yes") {
    saveRDS(data.scores, "./plots/dataFrames/plotData_plotScatter.rds")
  }

  # Return the plot. Adding more informative names describing plots
  out <- list("ellipseNMDS" = p1, "facetNMDS"= p2,"legend" = p3, "annualCentroids" = p4, 
              "annualTrajectories" =p5) # p3 is just the legend
  return(out)
}
