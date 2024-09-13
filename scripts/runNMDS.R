#####################################################################################################################
# Function to subset and/or aggregate full species accum matrix, and complete NMDS ordination
#####################################################################################################################

runNMDS <- function(accum, yrAgg, k, try, trymax, saveNMDS, preyCutoff, familyAgg, regionAgg) {
  accumSub <- accum

  # Now aggregate by year if required
  if(yrAgg == "yes" & regionAgg == "no") {
    accumYr <- aggregate(. ~ Year + Region + Predator_Species, accumSub[!names(accumSub) %in% c("Predator_ID")], 
                         FUN = mean)
    accumYr$count <- aggregate(. ~ Year + Region + Predator_Species, accumSub[!names(accumSub) %in% c("Predator_ID")], 
                               FUN = length)[,4]
    # Include at least N samples in a Region/Year
    accumYr <- subset(accumYr, count >= 1)
    accumYr$count <- NULL
    accumSub <- accumYr
  }
  if(yrAgg == "yes" & regionAgg == "yes") {
    accumYr <- aggregate(. ~ Year + Predator_Species, accumSub[!names(accumSub) %in% c("Predator_ID", "Region")], 
                         FUN = mean)
    accumYr$count <- aggregate(. ~ Year + Predator_Species, accumSub[!names(accumSub) %in% c("Predator_ID", "Region")], 
                               FUN = length)[,3]
    # Include at least N samples in a Region/Year
    accumYr <- subset(accumYr, count >= 1)
    accumYr$count <- NULL
    accumSub <- accumYr
  }

  #####################################################################################################################
  # Run NMDS. k is the no. dimensions to parse out and save
  nmds <- metaMDS(accumSub[!names(accumSub) %in% c("Predator_Species", "Predator_ID", "Year", "Region", "gutsPerSpPerYr")], 
                  distance = "horn", k = k, try = try, trymax = trymax, noshare = T)
  
  # Return/save the NMDS, and the subsetted/aggregated species matrix used to make it
  out <- list("nmds" = nmds, "accumSub" = accumSub)
  # Save object if specified. This is a list with both the NMDS and input data
  if(saveNMDS == "yes") {
    saveRDS(out, paste0("./data/nmds_", preyCutoff, "_", "familyAgg", familyAgg, "_yrAgg", yrAgg,
                        "_", min(accumSub$Year), "_", max(accumSub$Year), ".rds"))
  }
  return(out)
}