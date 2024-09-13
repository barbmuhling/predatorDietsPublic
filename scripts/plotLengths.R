#############################################################################################################
# Compare measured prey lengths across predators
# Note lengths are likely opportunistic, and may not be representative
#############################################################################################################

plotLengths <- function(hms_prey_lengths, paletteSpecies, plotMax, savePlot, saveDF) {
  # Only albacore, bluefin, swordfish, common thresher have > 100 prey measured
  # hms_prey_lengths_sub <- subset(hms_prey_lengths, Predator_Species == "Albacore" | Predator_Species == "Bluefin tuna" |
  #                                  Predator_Species == "Broadbill Swordfish" | Predator_Species == "Common Thresher Shark")
  # hms_prey_lengths_sub <- subset(hms_prey_lengths,)
  # Define which predators are present here
  pred_order <- c("Albacore", "Bluefin tuna", "Broadbill Swordfish",
                  "Blue Shark","Shortfin Mako Shark", "Common Thresher Shark", "Bigeye Thresher Shark",
                   "Short-Beaked Common Dolphin", "Long-Beaked Common Dolphin", "Northern Right Whale Dolphin")
  paletteSpecies$grp <- reorder.factor(paletteSpecies$grp, new.order=pred_order)
  paletteSpecies <-paletteSpecies[order(paletteSpecies$grp),]
  hms_prey_lengths$Predator_Species <- reorder.factor(hms_prey_lengths$Predator_Species, new.order=pred_order)
  
  cts <- aggregate(Prey_Length ~ Predator_Species, data = subset(hms_prey_lengths, Year < 2016), FUN = length)
  colnames(cts)[2] <- "meas_num"
  hms_prey_lengths_tot <- merge(hms_prey_lengths, cts)
  
  tots <- aggregate(Prey_N~Predator_Species, data = subset(hms_prey, Year<2016), FUN = sum)
  hms_prey_lengths_tot <- merge(hms_prey_lengths_tot, tots)
  hms_prey_lengths_tot$percent_measured <- round(hms_prey_lengths_tot$meas_num/hms_prey_lengths_tot$Prey_N, digits = 2)
  
  # make a list of labels to export so measurement count and percent data are readily available
  hms_prey_lengths_tot$labels <- as.factor(paste(hms_prey_lengths_tot$Predator_Species, as.character(hms_prey_lengths_tot$meas_num), 
                                                 as.character(hms_prey_lengths_tot$percent_measured), sep = "_"))
  meas2015 <- list(unique(hms_prey_lengths_tot$labels))
  
  # create ordered variable that shortens and reorders predator name
  hms_prey_lengths_tot$abs <- hms_prey_lengths_tot$Predator_Species
  levels(hms_prey_lengths_tot$abs) <- rev(list("ALBT"= "Albacore",
                                          "PBFT"= "Bluefin tuna",
                                          "NRWD"= "Northern Right Whale Dolphin",
                                          "LBCD"= "Long-Beaked Common Dolphin"  ,
                                          "SBCD"= "Short-Beaked Common Dolphin",
                                          "CMTS"= "Common Thresher Shark",
                                          "SFMS"= "Shortfin Mako Shark",
                                          "BLS"= "Blue Shark",
                                          "BBSF"= "Broadbill Swordfish" ,
                                          "BETS"= "Bigeye Thresher Shark"))
  # get color palette
  paletteSpecies2 <- subset(paletteSpecies, grp %in% hms_prey_lengths_tot$Predator_Species) 
  paletteSp <- paletteSpecies2$col # Groups predators by color 
  
  # Plot prey sizes for all predators 1998-2015
  hms_prey_lengths2015 <- subset(hms_prey_lengths_tot, Year<2016)
  p1 <- ggplot(hms_prey_lengths2015) + 
    geom_boxplot(aes(x = abs, y = Prey_Length, fill = Predator_Species, color = Predator_Species), alpha =0.65, size = 1.25) +
    scale_fill_manual("Predator Species", values = c(paletteSp))+ 
    scale_color_manual("Predator Species", values = c(paletteSp)) + theme_bw(base_size =16) +
    geom_hline(yintercept = 50, linetype = 2)+
    geom_hline(yintercept = 200, linetype = 2)+
    ylab("Prey Length (mm)")+
    coord_flip()+
    # facet_wrap(~group)+
    scale_y_continuous(limits = c(0, 875), expand = c(0.01,0))+
    theme( legend.position = "none",
          panel.grid = element_blank(),axis.title.y = element_blank())  
  p1
  
  # plot p1 without labels for figure assembly in Photoshop
  p1nolab<- p1+
    theme( legend.position = "none",axis.text=element_blank(),
           panel.grid = element_blank(),axis.title = element_blank())  
  p1nolab
  
  # summarize maximum prey size per stomach. Is perhaps more conservative given variability in prey measurements across predators
  hms_prey_lengths2015max <- aggregate(Prey_Length~Predator_Species*abs*Predator_ID*Prey_Length,hms_prey_lengths_tot, FUN = max)
  hms_pred_max_meas_count <- aggregate(Prey_Length~Predator_Species,hms_prey_lengths2015max, FUN = length )
  colnames(hms_pred_max_meas_count) <- c("Predator_Species","pred_count_wmeas")
  hms_pred_tot_count <- subset(unique(hms_prey[,c("Predator_ID", "Predator_Species")]))
  hms_pred_tot_count <-aggregate(Predator_ID~Predator_Species,hms_pred_tot_count, FUN = length)
  colnames(hms_pred_tot_count) <- c("Predator_Species","pred_count")
  
  pred_summary <- merge(hms_pred_max_meas_count , hms_pred_tot_count, by = "Predator_Species")
  pred_summary$percent_meas <- round(pred_summary$pred_count_wmeas/pred_summary$pred_count, 2)
    
  #make a list of labels to export so measurement count and percent data are readily available
  pred_summary$labels <- as.factor(paste(pred_summary$Predator_Species, as.character(pred_summary$pred_count_wmeas), as.character(pred_summary$percent_meas), sep = "_"))
  pred_measure_summ <- list(unique(pred_summary$labels))
  
  p1max <- ggplot(hms_prey_lengths2015max) + 
    geom_boxplot(aes(x = abs, y = Prey_Length, fill = Predator_Species, color = Predator_Species), alpha =0.65, size = 1.25) +
    scale_fill_manual("Predator Species", values = c(paletteSp)) + 
    scale_color_manual("Predator Species", values = c(paletteSp)) + theme_bw(base_size =16) +
    geom_hline(yintercept = 50, linetype = 2)+
    geom_hline(yintercept = 200, linetype = 2)+
    ylab("Max prey Length (mm)")+
    coord_flip()+
    scale_y_continuous(limits = c(0, 875), expand = c(0.01,0))+
    theme( legend.position = "none",
           panel.grid = element_blank(),axis.title.y = element_blank())  
  p1max
  
  # plot p1 without labels for figure assembly in Photoshop
  p1maxnolab<- p1max+
    theme( legend.position = "none",axis.text=element_blank(),
           panel.grid = element_blank(),axis.title = element_blank())  
  p1maxnolab
  # create a subset of just bluefin and swordfish
  sub_lengths2015 <- subset(hms_prey_lengths2015, (Predator_Species == "Bluefin tuna"|Predator_Species == "Broadbill Swordfish"))
  #assign a subset index for plotting temporal subset in correct order
  sub_lengths2015$subset <- "2.<2016"

  # create a 2016-018 subset for two target predators
  # need to recalculate total prey and percent measured for this new temporal subset
  sub_lengths2018 <- subset(hms_prey_lengths, (Predator_Species == "Bluefin tuna"|Predator_Species == "Broadbill Swordfish") &
                              Year > 2015)
  cts <- aggregate(Prey_Length~Predator_Species, data = sub_lengths2018, FUN = length)
  colnames(cts)[2] <- "meas_num"
  sub_lengths2018 <- merge(sub_lengths2018, cts)
  
  tots <- aggregate(Prey_N~Predator_Species, data = subset(hms_prey, Year > 2015), FUN = sum)
  sub_lengths2018 <- merge(sub_lengths2018, tots)
  sub_lengths2018$percent_measured <- round(sub_lengths2018$meas_num/sub_lengths2018$Prey_N, digits = 2)
  
  #make a list of labels to export so measurement count and percent data are readily available
  sub_lengths2018$labels <- as.factor(paste(sub_lengths2018$Predator_Species, as.character(sub_lengths2018$meas_num), 
                                            as.character(sub_lengths2018$percent_measured), sep = "_"))
  meas2018 <- list(unique(sub_lengths2018$labels))
  
  # create ordered variable that shortens and reorders predator name
    sub_lengths2018$abs <- sub_lengths2018$Predator_Species
  levels(sub_lengths2018$abs) <- rev(list("PBFT"= "Bluefin tuna",
                                             "BBSF"= "Broadbill Swordfish"))
  
  #assign a subset index for plotting temporal subset in correct order
  sub_lengths2018$subset <- "1.>=2016"
    
  #combine temporal and species subset dataframes together
  target_sub <- rbind(sub_lengths2015, sub_lengths2018)
  
  # split combined temporal subset dfs by predator to match formatting of Figure 6
  target_pbft <- subset(target_sub, abs == "PBFT")
  paletteSpecies3 <- subset(paletteSpecies, grp %in%  target_pbft$Predator_Species) 
  paletteSp <- paletteSpecies3$col # Groups predators by color 
  
  # plot variability in bluefin prey sizes across time periods
  p2pbft <- ggplot(target_pbft) + 
    geom_boxplot(aes(x = subset, y = Prey_Length, fill = Predator_Species,color = Predator_Species), alpha =0.65, size = 1.25) +
    scale_fill_manual("Predator Species", values = c(paletteSp)) + 
    scale_color_manual("Predator Species", values = c(paletteSp)) + theme_bw(base_size =16) +
    geom_hline(yintercept = 50, linetype = 2)+
    geom_hline(yintercept = 200, linetype = 2)+
    ylab("Prey Length (mm)")+
    coord_flip()+
    # facet_wrap(~group)+
    scale_y_continuous(limits = c(0, 665), expand = c(0.01,0))+
    theme( legend.position = "none",
           panel.grid = element_blank(),axis.title.y = element_blank())  
  p2pbft
  p2pbftNOlab <- p2pbft+
    theme( legend.position = "none",axis.text = element_blank(),
           panel.grid = element_blank(),axis.title= element_blank()) 
  p2pbftNOlab 
  
  # split combined temporal subset dfs by predator to match formatting of Figure 6
  target_bbsf <- subset(target_sub, abs == "BBSF")
  paletteSpecies3 <- subset(paletteSpecies, grp %in%  target_bbsf$Predator_Species) 
  paletteSp <- paletteSpecies3$col # Groups predators by color 
  
  # plot variability in swordfish prey sizes across time periods
    p2bbsf <- ggplot(target_bbsf) + 
    geom_boxplot(aes(x = subset, y = Prey_Length, fill = Predator_Species,color = Predator_Species), alpha =0.65, size = 1.25) +
    scale_fill_manual("Predator Species", values = c(paletteSp)) + 
    scale_color_manual("Predator Species", values = c(paletteSp)) + theme_bw(base_size =16) +
    geom_hline(yintercept = 50, linetype = 2)+
    geom_hline(yintercept = 200, linetype = 2)+
    ylab("Prey Length (mm)")+
    coord_flip()+
    # facet_wrap(~group)+
    scale_y_continuous(limits = c(0, 665), expand = c(0.01,0))+
    theme( legend.position = "none",
           panel.grid = element_blank(),axis.title.y = element_blank())  
    p2bbsf
    
    p2bbsfNOlab <- p2bbsf+
      theme( legend.position = "none",axis.text = element_blank(),
             panel.grid = element_blank(),axis.title= element_blank()) 
    p2bbsfNOlab 
  
  if(savePlot == "yes") {
    ggsave("./plots/preyLengthsAcrossPredators1998-2015.tiff", plot = p1, width = 4000, height = 2500, 
           units = "px", dpi = 400, compression = "lzw")
    ggsave("./plots/preyLengthsAcrossPredators1998-2015_NOlabs.tiff", plot = p1nolab, width = 4000, height = 2500, 
           units = "px", dpi = 400, compression = "lzw")
    ggsave("./plots/preyLengthsAcrossPredators1998-2015_MAXlengths.tiff", plot = p1max, width = 4000, height = 2500, 
           units = "px", dpi = 400, compression = "lzw")
    ggsave("./plots/preyLengthsAcrossPredators1998-2015_MAXlengths_NOlab.tiff", plot = p1maxnolab, width = 4000, height = 2500, 
           units = "px", dpi = 400, compression = "lzw")
    ggsave("./plots/preyLengths_Bluefin_time_comparison.tiff", plot = p2pbft, width = 4000, height = 500, 
           units = "px", dpi = 400, compression = "lzw")
    ggsave("./plots/preyLengths_Bluefin_time_comparison_NOlab.tiff", plot = p2pbftNOlab, width = 1500, height = 500, 
           units = "px", dpi = 400, compression = "lzw")
    ggsave("./plots/preyLengths_Swordfish_time_comparison.tiff", plot = p2bbsf, width = 4000, height = 500, 
           units = "px", dpi = 400, compression = "lzw")
    ggsave("./plots/preyLengths_Swordfish_time_comparison_NOlab.tiff", plot = p2bbsfNOlab, width = 1500, height = 500, 
           units = "px", dpi = 400, compression = "lzw")
  }
  if(saveDF == "yes") {
    saveRDS(hms_prey_lengths_sub, "./plots/dataFrames/plotData_plotLengths.rds")
  }
    out <- list("AllLengthsThrough2015" = p1, "MaxLengthsThrough2015" = p1max, "BluefinTemporal" = p2pbft, "SwordfishTemporal" = p2bbsf, 
                "MeasurementSummary2015"=meas2015, "MeasurementSummary2018"=meas2018, "MaxMeasSummary" = pred_measure_summ) # Return plots and objects all together in a list
    return(out)
}
