#################################################################################################################
# Create accumulation curves per species for full dataset 
# and curves per year, separately per species (one plot per species with curves per year included in annual NMDS)
#################################################################################################################

plotAccumulationCurves <- function(accum, paletteSpecies, type, savePlot, nboot, q) {
    # Set diet matrix to generate accumulation curves in iNEXT and subset
    # accum_curve <- accum
    accum_curve <- subset(accum, Year < 2016)
    
    # convert proportional abundance to occurrence
    accum_curve <- cbind(accum_curve[, 1:5], 
                         accum_curve[, 6:ncol(accum_curve)] %>% mutate_if(is.numeric, ~1 * (. != 0))) # Modification keeps Year
    
    #############################################################################################################
    # Calculate family richness per species per year, as well as total richness
    accumLong <- pivot_longer(accum_curve, cols = colnames(accum_curve[6:ncol(accum_curve)]),
                              names_to = "family", values_to = "presabs")
    accumAgg <- aggregate(presabs ~ Predator_Species + Year + family, accumLong, FUN = max, na.rm = TRUE)
    accumAggYr <- aggregate(presabs ~ Predator_Species + Year, accumAgg, FUN = sum, na.rm = TRUE)
    accumAggTotal <- aggregate(presabs ~ Predator_Species + family, accumLong, FUN = max, na.rm = TRUE)
    accumAggTotal <- aggregate(presabs ~ Predator_Species, accumAggTotal, FUN = sum, na.rm = TRUE)
    colnames(accumAggTotal)[2] <- "totalFamilies"
    accumAggYr <- left_join(accumAggYr, accumAggTotal, by = "Predator_Species")
    colnames(accumAggYr)[3] <- "familiesPerYear"
    # Annual richness as a proportion of total richness
    accumAggYr$annualPropnRichness <- accumAggYr$familiesPerYear / accumAggYr$totalFamilies
    # Add guts/sp/yr
    gutsSpYr <- aggregate(gutsPerSpPerYr ~ Predator_Species + Year, accum_curve, FUN = mean, na.rm = TRUE)
    accumAggYr <- left_join(accumAggYr, gutsSpYr, by = c("Predator_Species", "Year"))

    ###################################################################################################################
    # List of species to loop through (specifying instead of using unique to preserve desired order)
    # Also include the species abbreviations that are used later
    spNames <- data.frame("spName" = c("Albacore","Bluefin tuna", "Broadbill Swordfish", "Blue Shark", "Shortfin Mako Shark",
                 "Common Thresher Shark", "Bigeye Thresher Shark", "Short-Beaked Common Dolphin",
                 "Long-Beaked Common Dolphin", "Northern Right Whale Dolphin"),  
                          "spAbb" = c("alb","bft","swf","bls","sfm","cmt","bet","sbd","lbd","nwd"))
    
    # A function to process the data for each species so it can be plotted
    formatSpData <- function(spName) {
      # Note these col indices will only work if the leading columns added stay the same
      mySpData <- data.frame(subset(accum_curve, Predator_Species == spName)[,c(1,2,6:ncol(accum_curve))])
      # Change row names to be the specimen number (row names as sites is specification of iNEXT)
      rownames(mySpData) <- mySpData$Predator_ID
      # Remove column of specimen ids from matrix now that row names are the ids
      mySpData <- mySpData[,c(3:27)]
      # Transpose data frames (iNEXT needs columns as sites and rows as taxa)
      mySpData <- t(mySpData)
      return(mySpData)
    }
    
    # Loop through the species to generate 10 objects
    for(i in 1:nrow(spNames)) {
      spName <- spNames$spName[i]
      mySpData <- formatSpData(spName = spName)
      assign(x = spNames$spAbb[i], value = mySpData)
    }
    
    # Combine separate objects (of prey incidence per predator species) into a list
    # Reverse order to enable easier plotting below
    rarefication_test_groups <- rev(list(alb, bft, swf, bls, sfm, cmt, bet, sbd, lbd, nwd))
    
    # For iNEXT to treat different predators separately and run correctly, you must name objects in the list
    names(rarefication_test_groups) <- rev(spNames$spName) # Reverse order to enable easier plotting below
    
    # Run accumumlation curve function iNEXT datatype = incidence_raw using two different hill numbers, q = 0 and q =1
    # nboot = 500 for final analysis, can be lower for speed when testing
    if(q == 0) {
      interpolated_groups <- suppressWarnings(iNEXT(rarefication_test_groups, q = 0, datatype = "incidence_raw", 
                                                  nboot = nboot, endpoint = 1000)) 
      plotlabel <- "Species richness (Hill number of order q=0)"
    } else if (q == 1) {
      interpolated_groups <- suppressWarnings(iNEXT(rarefication_test_groups, q = 1, datatype = "incidence_raw", 
                                                   nboot = nboot, endpoint = 1000))
      plotlabel <- "Species evenness (Hill number of order q=1)"
    }
    
    # Run diversity estimation function
    diversity_estimates <- estimateD(rarefication_test_groups, datatype = "incidence_raw", base = "coverage",
                                     level = 0.95, conf = 0.95)
    
    # Save extrapolation values at endpoint for use in species-specific plots
    qds <- interpolated_groups$iNextEst$size_based
    qds <- subset(qds, t == 1000)
    # Add full sp names
    qds$spName <- qds$Assemblage
    qds <- left_join(qds, spNames, by = "spName")
    
    #######################################################################################################################
    ########################################### Make plots ################################################################
    #######################################################################################################################
    # Define which predators are present here
    paletteSpecies2 <- subset(paletteSpecies, grp %in% names(rarefication_test_groups))
    # Rename species column in palette to match working dataframe
    colnames(paletteSpecies2)[1] <- "Predator_Species"
    paletteSpecies2$Predator_Species <- factor(paletteSpecies2$Predator_Species, 
                              levels = (spNames$spName))
    # Rearrange to group predators by color
    paletteSpecies2 <- paletteSpecies2 %>% 
      arrange(factor(Predator_Species))
    paletteSp <- (paletteSpecies2$col) # Groups predators by color
      
    # This orders predators to get correct color assignments
    interpolated_groups$iNextEst$size_based$Assemblage <- factor(interpolated_groups$iNextEst$size_based$Assemblage, 
                                                                 levels = spNames$spName)
    # Plot and save if desired
    if(type == "total" | type == "both") { 
      # Fig S2. accum_curveulation curves iNEXT sample coverage only for specimens and taxa included in diet analyses
      # toggle between "type" = 1, 2, and 3 to see all sample size, diversity, and coverage plot combos from iNEXT package
      pAll <- ggiNEXT(interpolated_groups, type = 1, se = F) + theme_bw(base_size = 12)+
        theme(panel.grid.minor = element_blank(), legend.position = "none")+
        labs(x = "Stomachs with prey", y = plotlabel) +
        scale_color_manual(values = paletteSp )+
        scale_fill_manual(values = paletteSp )+
        scale_x_continuous(lim = c(-2, 700), expand = c(0,0))+
        scale_shape_manual(values = rep(19, 10), guide = "none")
      
      gb3 <- ggplot_build(pAll + theme(legend.text = element_text(size = 10), panel.grid = element_blank()))
      gb3$data[[1]]$size <- 6
      gb3$data[[1]]$stroke <- 1
      gb3$data[[2]]$linewidth <- 1
      gb3$data[[2]]$alpha <- 0.65
      gt3 <- ggplot_gtable(gb3)
      grid.draw(gt3)
      
      # pAll
      if(savePlot == "yes") {
         ggsave(paste0("./plots/accumCurves/Multipredator_iNEXT_accum_curveulation_q", q, "_noSE.jpg"),
                plot = gt3, width = 7, height = 4, units = "in", dpi = 400)
      }
  }
 
  #############################################################################################################
  ###################################### Plots per species per year ###########################################
  #############################################################################################################
  if(type == "bySpecies" | type == "both") {
      # Only for predators represented in at least 5 years of the time series by at least 3 individuals per year 
      # First define those predators
      preds <- aggregate(Predator_ID ~ Predator_Species + Year, accum, FUN = length)
      preds <- subset(preds, Predator_ID >= 3)
      preds2 <- aggregate(Predator_ID ~ Predator_Species, preds, FUN = length)
      preds2 <- subset(preds2, Predator_ID >= 5) # Leaves us with 8 predators
      
      # Add shorter names in, easier for list outputs
      preds2$spName <- preds2$Predator_Species
      preds2 <- left_join(preds2, spNames, by = "spName")
      
      # Build the iNEXT inputs, but this time by year within species
      accum_curve_sp <- accum
      
      # Setting year as factor, important for comparing within predators across years
      # Note: here we are not subsetting to 1998-2015
      accum_curve_sp$Year <- as.factor(accum_curve_sp$Year)
      
      # Convert proportional abundance to occurrence
      accum_curve_sp <- accum_curve_sp %>% mutate_if(is.numeric, ~1 * (. != 0))
    
      # A function to process the data for each species within year so it can be plotted
      formatSpDataYr <- function(spName, Yr) {
        # Note these col indices will only work if the leading columns added stay the same
        mySpDataYr <- data.frame(subset(accum_curve_sp, Predator_Species == spName & Year == Yr)[,c(1,2,6:ncol(accum_curve_sp))])
        # Remove prey taxa with all zeroes, likely causes issues later
        mySpDataYr2 <- mySpDataYr[, 3:ncol(mySpDataYr)] # Drop out the string cols just for colsums
        mySpDataYr2 <- mySpDataYr2[, colSums(mySpDataYr2) > 0] # Select prey with > 0 occurrences
        mySpDataYr <- cbind(mySpDataYr[, 1:2], mySpDataYr2) # Join back
        # Change row names to be the specimen number (row names as sites is specification of iNEXT)
        rownames(mySpDataYr) <- mySpDataYr$Predator_ID
        # Remove column of specimen ids from matrix now that row names are the ids
        mySpDataYr$Predator_Species <- mySpDataYr$Predator_ID <- mySpDataYr$Yr <- NULL
        # Transpose data frames (iNEXT needs columns as sites and rows as taxa)
        mySpDataYr <- t(mySpDataYr)
        return(mySpDataYr)
      }
      
      # Loop through the species and years to generate 8 lists
      for(j in 1:nrow(preds2)) {
        spName <- preds2$Predator_Species[j]
        # Use shorter names for lists or is too unwieldy
        spName2 <- preds2$spAbb[j]
        accum_curve_sub <- subset(accum_curve_sp, Predator_Species == spName)
        # Define which years have data available
        uniqueYrs <- aggregate(Predator_ID ~ Year, accum_curve_sub, FUN = length)
        uniqueYrs <- subset(uniqueYrs, Predator_ID >= 3) # Just years with >= 3 individuals 
        # Output this, helps with palette subsetting later
        assign(x = paste0(spName2, "UniqueYrs"), value = uniqueYrs)
        # Loop through years
        suppressWarnings(rm(rare, rareOut))
        for(k in 1:nrow(uniqueYrs)) {
           mySpDataYr <- formatSpDataYr(spName = spName, Yr = uniqueYrs$Year[k])
           # Don't include years with 5 or less prey taxa
           if(nrow(mySpDataYr) <= 5) {
             next
           }
           assign(x = paste0(spName2, uniqueYrs$Year[k]), value = mySpDataYr)
           rare <- eval(parse(text = paste0("list(", spName2, uniqueYrs$Year[k], " = mySpDataYr)"))) 
           if(k == 1) {
             rareOut <- rare
           } else {
             rareOut <- c(rareOut, rare) # We're combining matrices from each year within the loop
           }
        }
        assign(x = paste0(spName2, "_rare_test_groups_yr"), rareOut)
      }
      
      #######################################################################################################################
      ################################## Make species-specific plots ########################################################
      #######################################################################################################################
      # Add a function that can be applied to each predator
      savePlotsSpeciesYears <- function(rareTstGrp, spNameOut, spUniqueYrs) {
        # Endpoint is now calculated using the largest number of predators for any year
        # Means that endpoint is never more than twice the maximum no. samples per year
        endpoint <- round((max(sapply(rareTstGrp, ncol))) * 2, -1)
        
        # There is a bug/instability(?) in the iNEXT package that causes sbd to fail. The year causing the failure is 2014,
        # but it's not clear why. So we have to manually remove that year
        if(spNameOut == "sbd") {
          rareTstGrp$sbd2014 <- NULL
         }
        
        # This option typically works
        if(q == 0) {
          set.seed(123)
          interpolated_groupsYr <- suppressWarnings(iNEXT(rareTstGrp, q = 0, datatype = "incidence_raw", 
                                                        nboot = nboot, endpoint = endpoint))
          plotlabel <- "Species richness (Hill number of order q=0)"
          placeholderlimit <- 27
        } else if (q == 1) { # This one seems unstable
          set.seed(123)
          interpolated_groupsYr <- suppressWarnings(iNEXT(rareTstGrp, q = 1, datatype = "incidence_raw", 
                                                         nboot = nboot, endpoint = endpoint))
          plotlabel <- "Species evenness (Hill number of order q=1)"
          placeholderlimit <- 17
        }

        # Run diversity estimation function
        diversity_estimatesYr <- estimateD(rareTstGrp, datatype = "incidence_raw", base = "coverage",
                                           level = 0.95, conf = 0.95) 
        # interpolated_groupsYr$DataInfo
      
        # Define palette
        yrPalette <- data.frame("grp" = seq(1998, 2018), "col" = c("black", "black","black","black","black","black","black",
                                                                   "gray35", "gray35", "gray35", "gray35", "gray40", "gray40", 
                                                                   "gray40", "gray70", "gray70", "gray70", "gray70", "gray70",
                                                                   "gray70", "gray70"), 
                                "shape" = c(1,19,0,15,2,17,4,1,19,0,15,2,17,4,1,19,0,15,2,17,4))
        # Define which years are present here
        yrPalette2 <- subset(yrPalette, grp %in% spUniqueYrs$Year) 
        
        # Add color for hline
        paletteSpeciesSub <- paletteSpecies2 %>% 
          filter(Predator_Species==subset(qds, spAbb == spNameOut)$spName)
        paletteSp <- (paletteSpeciesSub$col) # Groups predators by color
        
        # Plot predator species curves by year. Add max qD from all-years analysis in as a ceiling
        # Note can use interpolated_groupsYr or interpolated_groups1Yr
        p2 <- ggiNEXT(interpolated_groupsYr, type = 1, se = F) + theme_bw(base_size = 12) + 
          geom_hline(data = subset(qds, spAbb == spNameOut), aes(yintercept = qD), linewidth = 4, alpha = 0.75, color = paletteSp) +
          theme(panel.grid.minor = element_blank(), legend.position = "right", legend.spacing.y = unit(0, "cm"))+
          labs(x = "Stomachs with prey", y = plotlabel) +
          scale_shape_manual(values = yrPalette2$shape,labels = spUniqueYrs$Year )+
          scale_color_manual(values = yrPalette2$col, labels = spUniqueYrs$Year )+
          scale_fill_manual(values = yrPalette2$col, labels = spUniqueYrs$Year )+
          # set flexible limits to include ideal amount of white space at x-axis minima for each plot
          scale_x_continuous(lim = c(-(endpoint/100), endpoint), expand = c(0,0))+
          scale_y_continuous(lim = c(0,placeholderlimit), expand = c(0,0))+
          theme(axis.title = element_blank(), 
                # axis.text = element_blank(), 
                legend.position = "none")
        
        # change size of points and width of lines
        gb3 <- ggplot_build(p2 + theme(legend.text = element_text(size = 10), panel.grid = element_blank()))
        gb3$data[[1]]$size <- 6
        gb3$data[[1]]$stroke <- 1
        gb3$data[[2]]$linewidth <- 1
        gb3$data[[2]]$alpha <- 0.65
        gt3 <- ggplot_gtable(gb3)
        grid.draw(gt3)
                  
        # Save
        if(savePlot == "yes") { 
          ggsave(paste0("./plots/accumCurves/", spNameOut, "Multipredator_iNEXT_accum_curveulation_q", q, "_noSE_blankAxes.jpg"), 
                 plot = gt3, width= 7, height = 4.5, units = "in", dpi = 400)
        }
        return(gt3)
      }
      
      # Go though predators to save plots
      # Often gives extrapolation warning, but given how we're calculating endpoint, this issue should be minimal in practice
      albp2 <- suppressWarnings(savePlotsSpeciesYears(rareTstGrp = alb_rare_test_groups_yr, spNameOut = "alb", 
                                     spUniqueYrs = albUniqueYrs))
      bftp2 <- suppressWarnings(savePlotsSpeciesYears(rareTstGrp = bft_rare_test_groups_yr, spNameOut = "bft", 
                                     spUniqueYrs = bftUniqueYrs))
      swfp2 <- suppressWarnings(savePlotsSpeciesYears(rareTstGrp = swf_rare_test_groups_yr, spNameOut = "swf", 
                                     spUniqueYrs = swfUniqueYrs))
      blsp2 <- suppressWarnings(savePlotsSpeciesYears(rareTstGrp = bls_rare_test_groups_yr, spNameOut = "bls", 
                                     spUniqueYrs = blsUniqueYrs))
      sfmp2 <- suppressWarnings(savePlotsSpeciesYears(rareTstGrp = sfm_rare_test_groups_yr, spNameOut = "sfm", 
                                     spUniqueYrs = sfmUniqueYrs))
      cmtp2 <- suppressWarnings(savePlotsSpeciesYears(rareTstGrp = cmt_rare_test_groups_yr, spNameOut = "cmt", 
                                     spUniqueYrs = cmtUniqueYrs))
      sbdp2 <- suppressWarnings(savePlotsSpeciesYears(rareTstGrp = sbd_rare_test_groups_yr, spNameOut = "sbd", 
                                     spUniqueYrs = sbdUniqueYrs))
      lbdp2 <- suppressWarnings(savePlotsSpeciesYears(rareTstGrp = lbd_rare_test_groups_yr, spNameOut = "lbd", 
                                     spUniqueYrs = lbdUniqueYrs))
  }
    
  # Return plots as list 
  if(type == "total") { # Return one plot with all species on it
    out <- list("allSp" = pAll)
  } else if (type == "bySpecies") { # Return all species plots
    out <- list("alb" = albp2, "bft" = bftp2, "swf" = swfp2, "bls" = blsp2, "sfm" = sfmp2, "cmt" = cmtp2, 
                "sbd" = sbdp2, "lbd" = lbdp2) 
  } else if (type == "both") {
    out <- list("allSp" = pAll, "alb" = albp2, "bft" = bftp2, "swf" = swfp2, "bls" = blsp2, "sfm" = sfmp2, "cmt" = cmtp2, 
                "sbd" = sbdp2, "lbd" = lbdp2) 
  }
 return(out)
}