#############################################################################################################
# Assign an inshore/offshore habitat to each prey taxa, based on catches from SWFSC CPS trawl surveys
# Return bar graph showing prey types eaten per predator and per year
#############################################################################################################

plotInshoreOffshoreHabitats <- function(preyCutoff, startYr, endYr, cropStart, cropEnd,spp, hms_prey, prey_taxonomy, 
                                        savePlot, saveDF) {
  # Generate a new accum matrix
  source("./scripts/getSpMatrix.R")
  accum <- getSpMatrix(hms_prey = hms_prey, prey_taxonomy = prey_taxonomy, startYr = startYr, endYr = endYr, endYrAc = cropEnd,
                           preyCutoff = preyCutoff, saveFile = "no", familyAgg = "no")
  accumLong <- pivot_longer(accum, cols = 6:ncol(accum), names_to = "Prey_LPT", values_to = "Propn")
  
  accumLong <- subset(accumLong, Predator_Species %in% spp)
  accumLong <- subset(accumLong, (Year > (cropStart-1) & Year < (cropEnd+1)))
  # Assign each prey as inshore or offshore from CPS survey trawl catches
  inOff <- aggregate(Prey_N ~ Prey_LPT + Prey_Class + inOff_habitat, hms_prey, FUN = length)
  inOff$Prey_N <- NULL

  # Join with prey matrix
  accumLong <- left_join(accumLong, inOff, by = c("Prey_LPT")) 
  # Print information on how many prey had inshore/offshore habitat information 
  accumLongPos <- subset(accumLong, Propn > 0)
  print(paste0("There are ", sum(!is.na(accumLongPos$inOff_habitat)), " of ",  nrow(accumLongPos), 
               " total prey records with inshore/offshore habitat information"))
  
  # Summarize by inOff_habitat within Predator_ID
  predSum <- aggregate(Propn ~ Predator_ID + inOff_habitat + Prey_Class, accumLong, FUN = sum, na.rm = TRUE)
  # Total propns within a predator usually > 0.75, but can be < 1 if many high-level ID prey in guts (e.g. Pisces, Decapodiformes)
  # So we assume that the inOff_habitat of IDed prey is proportional
  # Calculate the total proportion in each predator gut, and then standardize against that
  totalPropn <- aggregate(Propn ~ Predator_ID, predSum, FUN = sum, na.rm = TRUE)
  colnames(totalPropn)[2] <- "IDedPropn"
  predSum <- left_join(predSum, totalPropn, by = "Predator_ID")
  predSum$adjPropn <- predSum$Propn / predSum$IDedPropn
  
  # Add in Year and Region
  yrRegion <- aggregate(Year ~ Predator_ID + Region + gutsPerSpPerYr + Predator_Species, accumLong, FUN = mean, na.rm = TRUE)
  predSum <- left_join(predSum, yrRegion, by = c("Predator_ID"))
  
  # Now summarize across years
  predSumInOff <- aggregate(adjPropn ~ Predator_Species + Prey_Class + inOff_habitat, 
                            subset(predSum, gutsPerSpPerYr > 0), FUN = mean, na.rm = TRUE, na.action = na.pass) 
  # Make prey type name a little nicer
  predSumInOff$PreyType <- factor(paste0(predSumInOff$inOff_habitat, "_", predSumInOff$Prey_Class))
  levels(predSumInOff$PreyType) <- list("Inshore" = "inshore_Actinopterygii",
                                       "Offshore" = "offshore_Actinopterygii",
                                       "Inshore" = "inshore_Cephalopoda",
                                       "Offshore" = "offshore_Cephalopoda",
                                       "Inshore" = "inshore_Malacostraca",
                                       "Offshore" = "offshore_Thaliacea",
                                       "Inshore" = "inshore_Elasmobranchii",
                                       "Offshore" = "offshore_Elasmobranchii")
  # Re-order predator species
  predSumInOff$Predator_Species <- factor(predSumInOff$Predator_Species, 
                                          levels = rev(c("Common Thresher Shark","Long-Beaked Common Dolphin",
                                                         "Bluefin tuna", "Albacore",
                                                         "Blue Shark", "Broadbill Swordfish","Bigeye Thresher Shark", 
                                                         "Shortfin Mako Shark", "Short-Beaked Common Dolphin",
                                                         "Northern Right Whale Dolphin")))
  
  # Plot 
  p1 <- ggplot(predSumInOff) + 
    geom_bar(aes(x = Predator_Species, y = adjPropn, fill = PreyType), stat = "identity", position = "stack", alpha = 0.6) + 
    scale_fill_manual("Prey Type", values = c("#88D1B9", # Crustaceans
                                              "#2D2D9C")) + # Other
    ylab("Mean Proportion") + 
    theme_bw(base_size = 14) +
    coord_flip()+
    guides(fill = guide_legend(ncol = 1, byrow = TRUE))+
    theme(axis.title = element_blank(), panel.grid.major = element_blank(), legend.position = "bottom", 
          legend.title = element_blank(),panel.grid.minor = element_blank())
  p1
  # Plot without labels for easy assembly in Photoshop
  p2 <- ggplot(predSumInOff) + 
    geom_bar(aes(x = Predator_Species, y = adjPropn, fill = PreyType), stat = "identity", position = "stack", alpha = 0.6) + 
    scale_fill_manual("Prey Type", values = c("#88D1B9", # Crustaceans
                                              "#2D2D9C")) + # Other
    ylab("Mean Proportion") + 
    theme_bw(base_size = 30) +
    coord_flip(expand=0)+
    scale_y_continuous(position = "right")+
    theme(axis.title = element_blank(), axis.text = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), legend.position = "none", legend.title = element_blank())
  p2
  
  #get legend (requires ggpubr package)
  leg <- get_legend(p1)
  # Convert to a ggplot and print
  p3 <- as_ggplot(leg)
  
  if(savePlot == "yes") {
    ggsave(paste0("./plots/preyInshoreOffshoreBySpecies_SoCeCA_", cropStart, "_", cropEnd,  ".tiff"),
           plot = p1, width = 2500, height = 3000, units = "px", dpi = 400, compression = "lzw")
    ggsave(paste0("./plots/preyInshoreOffshoreBySpecies_SoCeCA_noLabel_", cropStart, "_", cropEnd,  ".tiff"),
           plot = p2, width = 1000, height = 3500, units = "px", dpi = 400, compression = "lzw")
    ggsave(paste0("./plots/preyInshoreOffshoreBySpecies_SoCeCA_LEGEND", cropStart, "_", cropEnd,  ".tiff"),
           plot = p3, width = 850, height = 1500, units = "px", dpi = 400, compression = "lzw")
  }
  if(saveDF == "yes") {
    saveRDS(predSumInOff, "./plots/dataFrames/plotData_plotInshoreOffshoreHabitats.rds")
  }
  out <- list("horizontal" = p1, "horizontalNoLabs" = p2, "legend" = p3) # Return plots all together in a list
  return(out)
}