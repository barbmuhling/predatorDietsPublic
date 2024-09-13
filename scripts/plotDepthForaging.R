#############################################################################################################
# Calculate the mean depth of foraging by day and night for each predator
# Can do by prey occurrence, or weight by prey abundance (Prey_N). Return bar graph
#############################################################################################################

plotDepthForaging <- function(accumDepth, hms_prey, preyMetric, hms_predator_depths, paletteSpecies, savePlot, saveDF) {
    # Drop out Predator_Species years with few guts sampled
    accumDepth <- subset(accumDepth, gutsPerSpPerYr > 0)
    # Pivot longer
    accumDepthLong <- pivot_longer(accumDepth, cols = 6:ncol(accumDepth), names_to = "Prey_LPT", values_to = "Propn")
    # Add presence/absence column
    accumDepthLong$pa <- ifelse(accumDepthLong$Propn > 0, 1, 0) # 78944
    # Add in depth_mode_day and depth_mode_night by predator-prey interaction
    accumDepthLong <- left_join(accumDepthLong, hms_prey[c("Predator_ID", "Prey_LPT", "depth_mode_day", "depth_mode_night")],
                                by = c("Predator_ID", "Prey_LPT")) 
    accumDepthLong <- accumDepthLong[!duplicated(accumDepthLong), ]
    # Use presence/absence, or abundance to estimate foraging depths by day and night
    if(preyMetric == "occurrence") {
      # Simple mean within predator_ID
      accumDepthLong$depth_mode_both <- (accumDepthLong$depth_mode_day+accumDepthLong$depth_mode_night)/2
       depthMean <- aggregate(cbind(depth_mode_day, depth_mode_night, depth_mode_both) ~ Predator_ID + Predator_Species, 
                               subset(accumDepthLong, pa == 1), FUN = mean, na.rm = TRUE)
    } else if (preyMetric == "abundance") {
        # Here we need the Prey_N 
        preyN <- aggregate(Prey_N ~ Predator_ID + Prey_LPT, hms_prey, FUN = mean, na.rm = TRUE)
        accumDepthLong <- left_join(accumDepthLong, preyN, by = c("Predator_ID", "Prey_LPT"))
        # Set NA to 0
        accumDepthLong$Prey_N <- ifelse(is.na(accumDepthLong$Prey_N), 0, accumDepthLong$Prey_N)
        # Calculate weighted mean depth
        accumDepthLong$Prey_Nday <- accumDepthLong$Prey_N * accumDepthLong$depth_mode_day
        accumDepthLong$Prey_Nnight <- accumDepthLong$Prey_N * accumDepthLong$depth_mode_night
        accumDepthLong$Prey_Nboth <- accumDepthLong$Prey_N * ((accumDepthLong$depth_mode_day+accumDepthLong$depth_mode_night)/2)
        depthMean <- aggregate(cbind(Prey_Nday, Prey_N) ~ Predator_ID + Predator_Species, accumDepthLong, FUN = sum, na.rm = TRUE)
        depthMean$depth_mode_day <- depthMean$Prey_Nday / depthMean$Prey_N
        # Add night
        nightMean <- aggregate(cbind(Prey_Nnight, Prey_N) ~ Predator_ID + Predator_Species, accumDepthLong, FUN = sum, na.rm = TRUE)
        nightMean$depth_mode_night <- nightMean$Prey_Nnight / nightMean$Prey_N
        depthMean <- left_join(depthMean, nightMean[c("Predator_ID", "Predator_Species","depth_mode_night")], 
                               by = c("Predator_ID", "Predator_Species"))
        # Add mean of day and night
        bothMean <- aggregate(cbind(Prey_Nboth, Prey_N) ~ Predator_ID + Predator_Species, accumDepthLong, FUN = sum, na.rm = TRUE)
        bothMean$depth_mode_both <- bothMean$Prey_Nboth / bothMean$Prey_N
        depthMean <- left_join(depthMean, bothMean[c("Predator_ID", "Predator_Species","depth_mode_both")], 
                               by = c("Predator_ID", "Predator_Species"))
        depthMean$Prey_Nday <- depthMean$Prey_N <- NULL
    } else {
      print("Error in preyMetric: must be \"occurrence\" or \"abundance\"")
      stop()
    }
 
    #####################################################################################################
    # Easier to pivot longer by day/night for plotting
    predDepthLong <- pivot_longer(depthMean, cols = c("depth_mode_day", "depth_mode_night", "depth_mode_both"), 
                                                  names_to = "dn", values_to = "preyDepth")
    # Add year
    predYr <- aggregate(Year ~ Predator_ID, hms_prey, FUN = max, na.rm = TRUE)
    predDepthLong <- left_join(predDepthLong, predYr, by = "Predator_ID")
    # Add predator max depths
    predDepthLong <- left_join(predDepthLong, hms_predator_depths[c("Predator_Species", "pred_max_depth_day", "pred_max_depth_night", 
                                                     "pred_depth_mode_day", "pred_depth_mode_night")], by = c("Predator_Species"))
    predDepthLong$pred_depth_max <- ifelse(predDepthLong$dn == "depth_mode_day", predDepthLong$pred_max_depth_day, 
                                           ifelse(predDepthLong$dn == "depth_mode_night",predDepthLong$pred_max_depth_night,
                                                  ifelse(predDepthLong$dn == "depth_mode_both" & 
                                                           (predDepthLong$pred_max_depth_day > predDepthLong$pred_max_depth_night), 
                                                            predDepthLong$pred_max_depth_day, predDepthLong$pred_max_depth_night)))
    # Add estimated foraging times
    predDepthLong <- merge(predDepthLong, hms_predator_depths[, c("Predator_Species", "diel_foraging")], all.x =T)
    # Add dummy variable to use as x-axis when making desired plots
    predDepthLong$base <- "base"
    # Subsetting by estimated foraging time
    predDepthLong2 <- predDepthLong %>%
                             filter((diel_foraging == "night" & dn=="depth_mode_night")|
                                      (diel_foraging == "day" & dn=="depth_mode_day")|
                                      (diel_foraging == "both" & dn=="depth_mode_both")|
                                      Predator_Species == "Northern Right Whale Dolphin" & dn == "depth_mode_night")
    predDepthLong2$pred_depth_max <- ifelse(predDepthLong2$diel_foraging == "both" & 
                                              predDepthLong2$pred_max_depth_day > predDepthLong2$pred_max_depth_night, 
                                            predDepthLong2$pred_max_depth_day,
                                           ifelse(predDepthLong2$diel_foraging == "both" & 
                                                    predDepthLong2$pred_max_depth_day < predDepthLong2$pred_max_depth_night, 
                                                  predDepthLong2$pred_max_depth_night, predDepthLong2$pred_depth_max))
    
    # Subsetting by estimated foraging time
    predDepthLong3 <- predDepthLong %>%
      filter((diel_foraging == "night" & dn == "depth_mode_night")|
               (diel_foraging == "day" & dn == "depth_mode_day")|
               (diel_foraging == "both" & (dn == "depth_mode_day"|dn == "depth_mode_night"))|
               (Predator_Species == "Northern Right Whale Dolphin"& dn == "depth_mode_night"))
    predDepthLong3$dm <- ifelse(predDepthLong3$dn == "depth_mode_day", "day", 
                                ifelse(predDepthLong3$dn == "depth_mode_night", "night", NA))
    
    #####################################################################################################
    # Plot: Define which predators are present here
    schoeners_order <- c("Common Thresher Shark","Long-Beaked Common Dolphin",
                         "Bluefin tuna", "Albacore",
                         "Blue Shark", "Broadbill Swordfish","Bigeye Thresher Shark", "Shortfin Mako Shark",
                         "Short-Beaked Common Dolphin","Northern Right Whale Dolphin")
    paletteSpecies$grp <- reorder.factor(paletteSpecies$grp, new.order = schoeners_order)
    paletteSpecies <-paletteSpecies[order(paletteSpecies$grp),]
    paletteSp <- paletteSpecies$col # Groups predators by color 
    
    predDepthLong$Predator_Species <- reorder.factor(predDepthLong$Predator_Species, new.order=schoeners_order)
    predDepthLong2$Predator_Species <- reorder.factor(predDepthLong2$Predator_Species, new.order=schoeners_order)
    predDepthLong3$Predator_Species <- reorder.factor(predDepthLong3$Predator_Species, new.order=schoeners_order)
    
    # Plot summarizing estimated foraging depths based on foraging time of day 
    # (if "both" is an average of day and night for estimated foraging depth and max tag depth)
    # one LBCD does not have prey depth and NRWD does not have tag data for hline. both result in plot warnings.
    p2 <- ggplot(predDepthLong2) + # Can try occurrence or abundance here 
      geom_violin(aes(x = base, y = -preyDepth, group = dn, fill = Predator_Species), alpha = 0.5, 
                  draw_quantiles = 0.5, size = 1.125) + 
      facet_wrap(~ Predator_Species , nrow=1) +
      scale_fill_manual("Predator Species", values = c(paletteSp)) + 
      scale_color_manual("Predator Species", values = c(paletteSp)) +
      geom_hline(aes(yintercept = -pred_depth_max, color =Predator_Species), linewidth = 4,linetype=1, alpha = 1) +
      theme_bw(base_size = 18) +
      scale_y_continuous(expand = c(0,0))+
      guides(fill = guide_legend(ncol = 2, byrow = TRUE))+
      theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            strip.background = element_blank(),strip.text.x = element_blank(), legend.title = element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") 
    p2
    
    p3 <- ggplot(predDepthLong2) + # Can compare occurrence or abundance here 
      geom_violin(aes(x = base, y = -preyDepth, group = dn, fill = Predator_Species), alpha = 0.5, draw_quantiles = 0.5) + 
      facet_wrap(~ Predator_Species+diel_foraging , nrow=1) +
      scale_fill_manual("Predator Species", values = c(paletteSp)) + 
      scale_color_manual("Predator Species", values = c(paletteSp)) +
      geom_hline(aes(yintercept = -pred_depth_max, color =Predator_Species), linewidth = 4,linetype=1, alpha = 0.85) +
      theme_bw(base_size = 18) +
      scale_y_continuous(expand = c(0,0))+
      guides(fill = guide_legend(ncol = 2, byrow = TRUE))+
      theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            legend.title = element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") 
    p3
    
    # Plot summarizing estimated foraging depths based on foraging time of day 
    # (if "both" shows both day and night separately for estimated foraging depth and max tag depth)
    p4 <- ggplot(predDepthLong3) + # Can try occurrence or abundance here 
      geom_violin(aes(x = base, y = -preyDepth, group = dn, fill = Predator_Species), alpha = 0.5, draw_quantiles = 0.5) + 
      facet_wrap(~ Predator_Species+ dm, nrow=1) +
      scale_fill_manual("Predator Species", values = c(paletteSp)) + 
      scale_color_manual("Predator Species", values = c(paletteSp)) +
      geom_hline(aes(yintercept = -pred_depth_max, color =Predator_Species), linewidth = 4,linetype=1, alpha = 0.85) +
      theme_bw(base_size = 18) +
      scale_y_continuous(expand = c(0,0))+
      guides(fill = guide_legend(ncol = 2, byrow = TRUE))+
      theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            strip.background = element_blank(),strip.text.x = element_blank(), 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")
    p4
    
    # Guide for foraging times for p3
    p5 <- ggplot(predDepthLong3) + # Can compare occurrence or abundance here 
      geom_violin(aes(x = base, y = -preyDepth, group = dn, fill = Predator_Species), alpha = 0.5, draw_quantiles = 0.5) + 
      facet_wrap(~ Predator_Species+dm , nrow=1) +
      scale_fill_manual("Predator Species", values = c(paletteSp)) + 
      scale_color_manual("Predator Species", values = c(paletteSp)) +
      geom_hline(aes(yintercept = -pred_depth_max, color =Predator_Species), linewidth = 4,linetype=1, alpha = 0.85) +
      theme_bw(base_size = 18) +
      scale_y_continuous(expand = c(0,0))+
      guides(fill = guide_legend(ncol = 2, byrow = TRUE))+
      theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") 
    p5
    
    if(savePlot == "yes") {
      ggsave("./plots/depthOfForaging_bothDiel_avgTagDepth.tiff", plot = p2, width = 5000, height = 2650, 
             units = "px", dpi = 400, compression = "lzw")
      ggsave("./plots/depthOfForaging_bothDiel_avgTagDepthKEY.tiff", plot = p3, width = 5000, height = 3000, 
             units = "px", dpi = 400, compression = "lzw")
      ggsave("./plots/depthOfForaging_separateDiel.tiff", plot = p4, width = 5600, height = 2400, 
             units = "px", dpi = 400, compression = "lzw")
      ggsave("./plots/depthOfForaging_separateDielKEY.tiff", plot = p5, width = 5000, height = 3000, 
             units = "px", dpi = 400, compression = "lzw")
    }
    if(saveDF == "yes") {
      saveRDS(predDepthLong, "./plots/dataFrames/plotData_plotDepthofForaging.rds")
    }
    return(p2)
}