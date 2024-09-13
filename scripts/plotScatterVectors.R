#################################################################################################################
# This function plots a scatterplot of an NMDS with species vectors
# Function requires an NMDS object, and a species matrix.
#################################################################################################################

plotScatterVectors <- function(myNMDS, myAccum, preyCutoff, familyAgg, yrAgg, paletteSpecies, 
                               savePlot, saveDF, hypCut) {
  # Using the scores function from vegan to extract the site scores and convert to a data.frame
  data.scores <- as.data.frame(scores(myNMDS) [[1]])  ### I added list index, required
  # Create a column of site names, from the rownames of data.scores
  data.scores$site <- rownames(data.scores) 
  # Add the grp variable created earlier
  data.scores$grp <- myAccum$Predator_Species
  # Add region and year too
  data.scores$region <- myAccum$Region
  data.scores$year <- myAccum$Year
  # head(data.scores)
  
  # Define which predators are present here
  paletteSpeciesSub <- subset(paletteSpecies, grp %in% data.scores$grp) 
  paletteSp <- paletteSpecies$col # Groups predators by color 
  
  # Adding labels to family groups in ordination vectors that contribute > X dissim
  # Get the species matrix without ID cols
  spMat <- as.matrix(myAccum[!names(myAccum) %in% c("Predator_Species", "Predator_ID", "Year", "Region", 
                                                    "gutsPerSpPerYr")]) # exclude ID columns)
  # Fit factors to existing NMDS ordination
  NMDS.fam.fit <- envfit(myNMDS, spMat, permutations = 999)
  # Extract scores
  NMDS.fam.scrs <- as.data.frame(scores(NMDS.fam.fit, display = "vectors")) # save species intrinsic values into dataframe
  NMDS.fam.scrs <- cbind(NMDS.fam.scrs, Prey_LPT = rownames(NMDS.fam.scrs)) # add species names to dataframe
  NMDS.fam.scrs <- cbind(NMDS.fam.scrs, pval = NMDS.fam.fit$vectors$pvals) #add pvalues to dataframe so can select sign. species
  NMDS.fam.scrs <- cbind(NMDS.fam.scrs, abrev = abbreviate(NMDS.fam.scrs$Prey_LPT, minlength = 3)) # abbreviate species names
  NMDS.fam.scrs <- subset(NMDS.fam.scrs, pval <= 0.05) # subset data to show species significant at e.g. 0.05
  NMDS.fam.scrs$hyp <- sqrt(NMDS.fam.scrs$NMDS1 ^ 2 + NMDS.fam.scrs$NMDS2 ^ 2)
  NMDS.fam.scrs2 <- subset(NMDS.fam.scrs, hyp >= hypCut) # Cutoff for important prey, imp on either axis
  
  # Scatterplot with vectors
  pScatterVector <- ggplot() + 
    geom_segment(data = NMDS.fam.scrs2, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2), arrow = arrow(length = unit(0.25, "cm")), 
                 colour = "grey10", lwd = 0.35) + # add vector arrows of significant species
    ggrepel::geom_label_repel(data = NMDS.fam.scrs2, aes(x = NMDS1, y = NMDS2, label = Prey_LPT), cex = 3, direction = "both",
                              segment.size = 0, alpha = 0.85, max.overlaps = 100, force = 1, force_pull = 50) +
    scale_colour_manual("Predator Species", values = c(paletteSp)) + 
    scale_x_reverse()+
    geom_hline(yintercept = 0, linetype = 5)+
    geom_vline(xintercept = 0, linetype = 5)+
    theme_bw(base_size = 14) +
    theme(axis.title = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
          # axis.text = element_blank(), axis.ticks = element_blank())
  pScatterVector
  
  pScatterVector2 <- ggplot() + 
    geom_hline(yintercept = 0, linetype = 5)+
    geom_vline(xintercept = 0, linetype = 5)+
    geom_segment(data = NMDS.fam.scrs2, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2), arrow = arrow(length = unit(0.25, "cm")), 
                 colour = "grey10", lwd = 1, alpha =0.65) + # add vector arrows of significant species
    scale_colour_manual("Predator Species", values = c(paletteSp)) + 
    scale_x_reverse(limits = c(0.75,-0.75), breaks = c(0.5, 0, -0.5))+
    scale_y_continuous(limits = c(-0.45,0.425), breaks = c(-0.25, 0, 0.25))+
    theme_bw(base_size = 14) +
    theme(axis.title = element_blank(), panel.grid.major = element_blank(), 
          axis.text = element_blank(),
          panel.grid.minor = element_blank())
  pScatterVector2
  # Save
  if(savePlot == "yes") {
    ggsave(paste0("./plots/ordinations/scatterVectors_", preyCutoff, "_", "familyAgg", familyAgg, 
                  "_yrAgg", yrAgg, "_", min(myAccum$Year), "_", max(myAccum$Year), ".tiff"), 
           plot = pScatterVector, width = 2000, height = 1700, units = "px", dpi = 400, compression = "lzw")
    ggsave(paste0("./plots/ordinations/scatterVectorsNoLabel_", preyCutoff, "_", "familyAgg", familyAgg, 
                  "_yrAgg", yrAgg, "_", min(myAccum$Year), "_", max(myAccum$Year), ".tiff"), 
           plot = pScatterVector2, width = 2000, height = 1700, units = "px", dpi = 400, compression = "lzw")
  }
  if(saveDF == "yes") {
    saveRDS(list("data.scores" = data.scores, "NMDS.fam.scrs" = NMDS.fam.scrs), "./plots/dataFrames/plotData_plotScatterVectors.rds")
  }
  
  # Return the plots
  out <- list("vectorLabels" = pScatterVector, "vectors" = pScatterVector2)
  return(out)
}
