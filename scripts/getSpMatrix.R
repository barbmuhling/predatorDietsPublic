#################################################################################################################
# This function prepares diet data from multiple predators ready for multivariate analyses
# Function requires start and end years to use for trimming data, as well as cutoffs for rare species,
# and rare predators
#################################################################################################################

getSpMatrix <- function(hms_prey, prey_taxonomy, startYr, endYr, endYrAc, preyCutoff, saveFile, familyAgg, printOutputs = "no") {
  # Select years to include
  hms_prey_crop <- subset(hms_prey, Year >= startYr & Year <= endYr) 
  
  # If not aggregating to family:
  if(familyAgg == "no") {
    # Calculate total no. of prey per species in each predator gut, also fixes duplicate records of the same prey in one stomach
    hms_prey_sum <- aggregate(Prey_N ~ Predator_ID + Predator_Species + Prey_LPT, hms_prey_crop, FUN = sum) 
    
    # Calculate total number of prey in gut of any taxa (will be used for calculating proportional abundance)
    hms_prey_tot <- aggregate(Prey_N ~ Predator_ID, hms_prey_crop, FUN = sum) 
    colnames(hms_prey_tot)[2] <- "total_N"
    
    # Add count field (1) to get total number of guts with 1+ prey per predator species
    hms_prey_crop$count <- 1
    hms_pred_tot <- aggregate(count ~ Predator_Species, 
                              unique(hms_prey_crop[,c("count", "Predator_ID", "Predator_Species")]), FUN = sum)
    
    # Join field showing total no. prey per stomach to total prey data
    hms_prey_sum <- left_join(hms_prey_sum, hms_prey_tot, by = "Predator_ID")
    # Add field showing proportional abundance per stomach
    hms_prey_sum$Prop_N <- hms_prey_sum$Prey_N / hms_prey_sum$total_N
    # Calculate overall proportional contribution of each prey species to each predator 
    prey_species <- aggregate(Prop_N ~ Prey_LPT * Predator_Species, hms_prey_sum, FUN = sum)
    # Add the total number of predators with prey in guts per spp ( 1 value per predator sp)
    prey_species  <- left_join(prey_species, hms_pred_tot, by = "Predator_Species")
    # Calculate mean proportional prey abundance per predator species (e.g. anchovy is 10% of albacore diet on average)
    prey_species$mean_propn <- prey_species$Prop_N / prey_species$count
    
    #############################################################################################################
    ####################### Choose thresholds to subset out rare predator and prey species ######################
    #############################################################################################################
    # Only include prey species that contribute > X% in at least 1 predator
    prey_species_crop <- subset(prey_species, mean_propn >= preyCutoff)
    hms_pred_tot_crop <- hms_pred_tot
    # Remove prey species not identified to at least family
    preyClassOrder <- subset(prey_taxonomy, is.na(Prey_Family))[c("Prey_LPT")] 
    prey_species_crop <- subset(prey_species_crop, !Prey_LPT %in% preyClassOrder$Prey_LPT)
    
    # Subset to remove rare prey taxa and rare predators
    hms_prey_sum_cutoff <- subset(hms_prey_sum, Prey_LPT %in% prey_species_crop$Prey_LPT) 
    hms_prey_sum_cutoff <- subset(hms_prey_sum_cutoff, Predator_Species %in% hms_pred_tot_crop$Predator_Species)
  }
  
  # If aggregating to family:
  if(familyAgg == "yes") {
    # Calculate total no. of prey per species in each predator gut, also fixes duplicate records of the same prey in one stomach
    hms_prey_sum <- aggregate(Prey_N ~ Predator_ID + Predator_Species + Prey_Family, hms_prey_crop, FUN = sum) 
    
    # Calculate total number of prey in gut of any taxa (will be used for calculating proportional abundance)
    hms_prey_tot <- aggregate(Prey_N ~ Predator_ID, hms_prey_crop, FUN = sum) 
    colnames(hms_prey_tot)[2] <- "total_N"
    
    # Add count field (1) to get total number of guts with 1+ prey per predator species
    hms_prey_crop$count <- 1
    hms_pred_tot <- aggregate(count ~ Predator_Species, 
                              unique(hms_prey_crop[,c("count", "Predator_ID", "Predator_Species")]), FUN = sum)
    
    # Join field showing total no. prey per stomach to total prey data
    hms_prey_sum <- left_join(hms_prey_sum, hms_prey_tot, by = "Predator_ID")
    # Add field showing proportional abundance per stomach
    hms_prey_sum$Prop_N <- hms_prey_sum$Prey_N / hms_prey_sum$total_N
    # Calculate overall proportional contribution of each prey species to each predator 
    prey_species <- aggregate(Prop_N ~ Prey_Family * Predator_Species, hms_prey_sum, FUN = sum)
    # Add the total number of predators with prey in guts per spp ( 1 value per predator sp)
    prey_species  <- left_join(prey_species, hms_pred_tot, by = "Predator_Species")
    # Calculate mean proportional prey abundance per predator species (e.g. anchovy is 10% of albacore diet on average)
    prey_species$mean_propn <- prey_species$Prop_N / prey_species$count
    
    ###############################################################################################################################
    ####################### Choose thresholds to subset out rare predator and prey species ########################################
    ###############################################################################################################################
    # Only include prey species that contribute > X% in at least 1 predator
    prey_species_crop <- subset(prey_species, mean_propn >= preyCutoff)
    hms_pred_tot_crop <- hms_pred_tot
    # Remove prey species not identified to at least family
    preyClassOrder <- subset(prey_taxonomy, is.na(Prey_Family))[c("Prey_LPT")] 
    prey_species_crop <- subset(prey_species_crop, !Prey_Family %in% preyClassOrder$Prey_Family)
    
    # Subset to remove rare prey taxa and rare predators
    hms_prey_sum_cutoff <- subset(hms_prey_sum, Prey_Family %in% prey_species_crop$Prey_Family) 
    hms_prey_sum_cutoff <- subset(hms_prey_sum_cutoff, Predator_Species %in% hms_pred_tot_crop$Predator_Species)
  }
  
  # Add year and region 
  yrs <- aggregate(Year ~ Predator_ID + Region, hms_prey, FUN = mean, na.rm = TRUE)
  hms_prey_sum_cutoff <- left_join(hms_prey_sum_cutoff, yrs, by = c("Predator_ID"))

  # Calculate stomachs per species per year
  predators <- aggregate(Prey_N ~ Predator_Species + Predator_ID + Year, hms_prey_sum_cutoff, FUN = length)
  yrs2 <- aggregate(Prey_N ~ Predator_Species + Year, predators, FUN = length)
  colnames(yrs2)[3] <- "gutsPerSpPerYr" 
  hms_prey_sum_cutoff <- left_join(hms_prey_sum_cutoff, yrs2, by = c("Predator_Species", "Year"))
  
  # If making a temporal subset on full dataset
  if(endYrAc != endYr) {
    hms_prey_sum_cutoff<- subset(hms_prey_sum_cutoff, Year >= startYr & Year <= endYrAc)
  }
  # Pivot to wide format to show prey matrix, values are proportional abundance per gut
  # By species
  if(familyAgg == "no") {
  accum <- pivot_wider(subset(hms_prey_sum_cutoff, gutsPerSpPerYr > 0), 
                       id_cols = c(Predator_Species, Predator_ID, Year, Region, gutsPerSpPerYr), 
                       names_from = Prey_LPT, values_from = Prop_N)
  }
  # Or by family
  if(familyAgg == "yes") {
    accum <- pivot_wider(subset(hms_prey_sum_cutoff, gutsPerSpPerYr > 0), 
                         id_cols = c(Predator_Species, Predator_ID, Year, Region, gutsPerSpPerYr), 
                         names_from = Prey_Family, values_from = Prop_N)
  }
  
  # Fill 0s for NAs
  accum[is.na(accum)] <- 0
    
  # Export, filename includes cutoffs
  if(saveFile == "yes") {
    saveRDS(accum, file = paste0("./data/spMatrix_", preyCutoff, "_", "familyAgg", familyAgg, 
                                 "_", startYr, "_", endYr, ".rds"))
    print("savedFile")
  }
  
  # If desired, print some simple info about the spp matrix: first no. taxa 
  if(printOutputs == "yes") {
    nms <- data.frame("col" = colnames(accum))
    nTaxa <- nrow(subset(nms, col != "Predator_Species" & col != "Predator_ID" & col != "Year" & col != "Region" & 
                           col != "gutsPerSpPerYr"))
    print(paste0("There are ", nTaxa, " prey taxa in matrix"))
    # Then total numbers of prey
    print(paste0("There are a total of ", sum(hms_prey_sum_cutoff$Prey_N, na.rm = TRUE), " prey items in matrix"))
  }
  # Return matrix
  return(accum)
}