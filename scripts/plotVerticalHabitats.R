#############################################################################################################
# Assign a vertical habitat to each prey taxa, based on the size eaten by each predator
# Return bar graph showing prey types eaten per predator and per year
#############################################################################################################

plotVerticalHabitats <- function(preyCutoff, startYr, endYr,  cropStart, cropEnd,spp, hms_prey, prey_taxonomy, savePlot, saveDF) {
  # Generate a new accum matrix based on start and end years 
  source("./scripts/getSpMatrix.R")
  accum <- getSpMatrix(hms_prey = hms_prey, prey_taxonomy = prey_taxonomy, startYr = startYr, endYr = endYr, endYrAc = cropEnd,
                           preyCutoff = preyCutoff, saveFile = "no", familyAgg = "no")
  accumLong <- pivot_longer(accum, cols = 6:ncol(accum), names_to = "Prey_LPT", values_to = "Propn")
  
  accumLong <- subset(accumLong, Predator_Species %in% spp)
  accumLong <- subset(accumLong, (Year > (cropStart-1) & Year < (cropEnd+1)))
  # Then need to add the likely lifestage that each prey was consumed at by each predator
  # (not all prey were measured, so it's just one value for each predator-prey combination)
  lifestages <- aggregate(Prey_N ~ Predator_Species + Prey_LPT + Prey_Class + life_stage_final + vert_habitat, 
                          hms_prey, FUN = length) 
  # Note: this will drop any prey record where vert_habitat is NA, such as for high-level IDs (e.g. Pisces, Decapodiformes)
  lifestages$Prey_N <- NULL
  # Join with prey matrix
  accumLong <- left_join(accumLong, lifestages, by = c("Prey_LPT", "Predator_Species")) 
  # Print information on how many prey had vertical habitat information
  accumLongPos <- subset(accumLong, Propn > 0)
  print(paste0("There are ", sum(!is.na(accumLongPos$vert_habitat)), " of ",  nrow(accumLongPos), 
                              " total prey records with vertical habitat information"))
  
  # Summarize by vert_habitat within Predator_ID
  predSum <- aggregate(Propn ~ Predator_ID + vert_habitat + Prey_Class, accumLong, FUN = sum, na.rm = TRUE)
  # Total propns within a predator are usually > 0.75, but can be < 1 if high-level ID prey in guts (e.g. Pisces, Decapodiformes)
  # So we assume that the vert_habitat of IDed prey is proportional
  # Calculate the total proportion in each predator gut, and then standardize against that
  totalPropn <- aggregate(Propn ~ Predator_ID, predSum, FUN = sum, na.rm = TRUE)
  colnames(totalPropn)[2] <- "IDedPropn"
  predSum <- left_join(predSum, totalPropn, by = "Predator_ID")
  predSum$adjPropn <- predSum$Propn / predSum$IDedPropn
  
  # Add in Year and Region
  yrRegion <- aggregate(Year ~ Predator_ID + Region + gutsPerSpPerYr + Predator_Species, accumLong, FUN = mean, na.rm = TRUE)
  predSum <- left_join(predSum, yrRegion, by = c("Predator_ID"))
  
  # Now summarize across years
  predSumVert <- aggregate(adjPropn ~ Predator_Species + Prey_Class + vert_habitat, 
                           subset(predSum, gutsPerSpPerYr > 0), FUN = mean, 
                           na.rm = TRUE, na.action = na.pass) 
  # Make prey type name a little nicer
  predSumVert$PreyType <- factor(paste0(predSumVert$vert_habitat, "_", predSumVert$Prey_Class))
  levels(predSumVert$PreyType) <- list("Epipelagic" = "epipelagic_Actinopterygii",
                                       "Mesopelagic" = "mesopelagic_Actinopterygii",
                                       "Demersal"  = "demersal_Actinopterygii",
                                       "Epipelagic" = "epipelagic_Cephalopoda",
                                       "Mesopelagic" = "mesopelagic_Cephalopoda",
                                       "Epipelagic" = "epipelagic_Malacostraca",
                                       "Demersal" = "demersal_Cephalopoda",
                                       "Mesopelagic" = "mesopelagic_Malacostraca",
                                       "Demersal" = "demersal_Malacostraca",
                                       "Demersal" = "benthic_Actinopterygii",
                                       "Demersal" = "benthic_Cephalopoda",
                                       "Epipelagic" = "epipelagic_Mammalia",
                                       "Epipelagic" = "epipelagic_Thaliacea",
                                       "Demersal" = "demersal_Elasmobranchii",
                                       "Epipelagic" = "epipelagic_Elasmobranchii")
  # Re-order predator species
  predSumVert$Predator_Species <- factor(predSumVert$Predator_Species, 
                                      levels = rev(c("Common Thresher Shark","Long-Beaked Common Dolphin",
                                                     "Bluefin tuna", "Albacore",
                                                     "Blue Shark", "Broadbill Swordfish","Bigeye Thresher Shark", 
                                                     "Shortfin Mako Shark", "Short-Beaked Common Dolphin",
                                                     "Northern Right Whale Dolphin")))
  
  # Plot 
  p1 <- ggplot(predSumVert) + 
    geom_bar(aes(x = Predator_Species, y = adjPropn, fill = PreyType), stat = "identity", position = "stack", alpha = 0.6) + 
    scale_fill_manual("Prey Type", values = c("#D9D9D9", "#8C8C8C","#404040")) + # Other
    ylab("Mean Proportion") + 
    theme_bw(base_size = 14) +
    coord_flip()+
    guides(fill = guide_legend(ncol = 1, byrow = TRUE))+
    scale_y_reverse()+
    theme(axis.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          legend.position = "bottom", legend.title = element_blank())
  
  # Plot without labels for easy assembly in Photoshop
  p2 <- ggplot(predSumVert) + 
    geom_bar(aes(x = Predator_Species, y = adjPropn, fill = PreyType), stat = "identity", position = "stack", alpha = 0.6) + 
    scale_fill_manual("Prey Type", values = c("#D9D9D9", "#8C8C8C","#404040")) + # Other
    ylab("Mean Proportion") + 
    theme_bw(base_size = 30) +
    coord_flip(expand = 0)+
    scale_y_reverse()+
    scale_y_continuous(position = "right")+
    guides(fill = guide_legend(nrow = 2))+
    theme(axis.title = element_blank(), axis.text = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), legend.position = "none", legend.title = element_blank())
  p2
  
  # Get legend (requires ggpubr package)
  leg <- get_legend(p1)
  # Convert to a ggplot and print
  p3 <- as_ggplot(leg)
  # Save
    if(savePlot == "yes") {
      ggsave(paste0("./plots/preyVerticalHabitat_SoCeCA_", cropStart, "_", cropEnd,  ".tiff"),
             plot = p1, width = 2500, height = 3000, units = "px", dpi = 400, compression = "lzw")
      ggsave(paste0("./plots/preyVerticalHabitat_SoCeCA_noLabel_", cropStart, "_", cropEnd,  ".tiff"),
             plot = p2, width = 1000, height = 3500, units = "px", dpi = 400, compression = "lzw")
      ggsave(paste0("./plots/preyVerticalHabitat_SoCeCA_LEGEND", cropStart, "_", cropEnd,  ".tiff"),
             plot = p3, width = 850, height = 1500, units = "px", dpi = 400, compression = "lzw")
    }
    if(saveDF == "yes") {
      saveRDS(predSumInOff, "./plots/dataFrames/plotData_plotInshoreOffshoreHabitats.rds")
    }
    out <- list("horizontal" = p1, "horizontalNoLabs" = p2, "legend" = p3) # Return plots all together in a list
    return(out)
  }