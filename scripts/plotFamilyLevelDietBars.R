#############################################################################################################
# Return bar graph showing prey families eaten per predator
#############################################################################################################

plotDietBars <- function(preyCutoff, startYr, endYr, cropStart, cropEnd, hms_prey, prey_taxonomy,spp, savePlot, saveDF) {
  # Generate a new species  matrix. Only differs by start and end years of data
  source("./scripts/getSpMatrix.R")
  accum <- getSpMatrix(hms_prey = hms_prey, prey_taxonomy = prey_taxonomy, startYr = startYr, endYr = endYr, endYrAc = cropEnd,
                           preyCutoff = preyCutoff, saveFile = "no", familyAgg = "yes")
  accumLong <- pivot_longer(accum, cols = 6:ncol(accum), names_to = "Prey_LPT", values_to = "Propn")
  
  # Select specified predators
  accumLong <- subset(accumLong, Predator_Species %in% spp)
  accumLong <- subset(accumLong, (Year > (cropStart-1) & Year < (cropEnd+1)))
  # Summarize by vert_habitat within Predator_ID
  predSum <- aggregate(Propn ~ Predator_ID + Prey_LPT, accumLong, FUN = sum, na.rm = TRUE)
  # Total propns within a predator are usually > 0.75, but can be < 1 if high-level ID prey in guts (e.g. Pisces, Decapodiformes)
  # Calculate the total proportion in each predator gut, and then standardize against that
  totalPropn <- aggregate(Propn ~ Predator_ID, predSum, FUN = sum, na.rm = TRUE)
  colnames(totalPropn)[2] <- "IDedPropn"
  predSum <- left_join(predSum, totalPropn, by = "Predator_ID")
  predSum$adjPropn <- predSum$Propn / predSum$IDedPropn
  
  # Add in Year and Region
  yrRegion <- aggregate(Year ~ Predator_ID + Region + gutsPerSpPerYr + Predator_Species, accumLong, FUN = mean, na.rm = TRUE)
  predSum <- left_join(predSum, yrRegion, by = c("Predator_ID"))
  
  # Now summarize across years
  predSumFam <- aggregate(adjPropn ~ Predator_Species + Prey_LPT, subset(predSum, gutsPerSpPerYr > 0), FUN = mean, 
                           na.rm = TRUE, na.action = na.pass) 
  
  # Re-order predator species
  predSumFam$Predator_Species <- factor(predSumFam$Predator_Species, 
                                         levels = rev(c("Common Thresher Shark","Long-Beaked Common Dolphin",
                                                        "Bluefin tuna", "Albacore",
                                                        "Blue Shark", "Broadbill Swordfish","Bigeye Thresher Shark", 
                                                        "Shortfin Mako Shark", "Short-Beaked Common Dolphin",
                                                        "Northern Right Whale Dolphin")))
  
  predSumFam$Prey_LPT <- factor(predSumFam$Prey_LPT,
                                 levels = c("Bathylagidae","Carangidae","Clupeidae","Engraulidae","Merlucciidae",
                                            "Microstomatidae","Mugilidae" ,"Myctophidae","Paralepididae","Scomberesocidae",
                                            "Scombridae","Scopelarchidae","Sebastidae", "Synodontidae","Trachipteridae",   
                                            "Alloposidae","Argonautidae","Enoploteuthidae","Gonatidae","Histioteuthidae",
                                            "Loliginidae","Octopodidae" ,"Octopoteuthidae","Ommastrephidae","Onychoteuthidae",
                                            "Munididae"))
 
  predSumFam <- subset(predSumFam, adjPropn>0)
  # Define which predators are present here
  palettePrey2 <- subset(palettePrey, Prey_LPT %in% predSumFam$Prey_LPT) 
  # Generate an ordered color palette
  palettePr <- palettePrey2$col # Groups predators by color 
  
  # Plot proportional barplots with legend
  p1 <- ggplot(predSumFam) + 
    geom_bar(aes(x = Predator_Species, y = adjPropn, fill = Prey_LPT), stat = "identity", position = "stack", alpha = 0.8) + 
    scale_fill_manual("Prey Type", values = palettePr) + # Other
    ylab("Mean Proportion") + 
    theme_bw(base_size = 20) +
    coord_flip(expand = 0)+
    scale_y_reverse()+
    guides(fill = guide_legend(ncol = 5, byrow = TRUE))+
    theme(axis.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          legend.position = "bottom", legend.title = element_blank())
  p1
  
  # version without labels for easy assembly in Photoshop
  p2 <- ggplot(predSumFam) + 
    geom_bar(aes(x = Predator_Species, y = adjPropn, fill = Prey_LPT), stat = "identity", position = "stack", alpha = 0.8) + 
    scale_fill_manual("Prey Type", values = palettePr) + # Other
    ylab("Mean Proportion") + 
    theme_bw(base_size = 30) +
    coord_flip(expand = 0)+
    scale_y_reverse(position = "right")+
    theme(axis.title = element_blank(), axis.text = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", 
          legend.title = element_blank())
  p2
  
  # Get legend (requires ggpubr package)
  leg <- get_legend(p1)
  # Convert to a ggplot and print
  p3 <- as_ggplot(leg)
  
  if(savePlot == "yes") {
    ggsave(paste0("./plots/preyCompositionByFamily_SoCeCA_", cropStart, "_", cropEnd,  ".tiff"), 
           plot = p1, width = 4500, height = 3000, units = "px", dpi = 400, compression = "lzw")
    ggsave(paste0("./plots/preyCompositionByFamily_SoCeCA_NoLabels_", cropStart, "_", cropEnd,  ".tiff"),
           plot = p2, width = 5500, height = 3500, units = "px", dpi = 400, compression = "lzw")
    ggsave(paste0("./plots/preyCompositionByFamily_SoCeCA_LEGEND_", cropStart, "_", cropEnd,  ".tiff"),
           plot = p3, width = 5000, height = 1500, units = "px", dpi = 400, compression = "lzw")
  }
  if(saveDF == "yes") {
    saveRDS(predSumVert, "./plots/dataFrames/plotData_plotFamilyBars.rds")
  }
  out <- list("MeanDiet" = p1, "MeanDietNoLabs" = p2, "legend" = p3) # Return plots all together in a list
  return(out)
}  