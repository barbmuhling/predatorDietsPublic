#############################################################################################################
# Compare measured lengths across predators. EFL for swordfish, TL for mammals, FL for other species
# Note lengths are likely derived from operculum length for most tunas
#############################################################################################################

plotLengthsPredators <- function(hms_predator_lengths, paletteSpecies, savePlot, saveDF) {
  # Define which predators are present here
  paletteSpecies2 <- subset(paletteSpecies, grp %in% hms_predator_lengths$Predator_Species)
  # Rename species column in palette to match working dataframe
  colnames(paletteSpecies2)[1] <- "Predator_Species"
  # Assign colors and reorder length dataframe by median predator length
  paletteSpecies2$Predator_Species <- factor(paletteSpecies2$Predator_Species, 
                                             levels = rev(c("Albacore", "Bluefin tuna", "Broadbill Swordfish",
                                      "Blue Shark","Shortfin Mako Shark", "Common Thresher Shark", "Bigeye Thresher Shark",
                                      "Short-Beaked Common Dolphin", "Long-Beaked Common Dolphin", "Northern Right Whale Dolphin")))
  paletteSpecies2 <- paletteSpecies2 %>% 
    arrange(factor(Predator_Species))
  paletteSp <- (paletteSpecies2$col)# Groups predators by color
  # Order by predator median lengths
  hms_predator_lengths$Predator_Species <- factor(hms_predator_lengths$Predator_Species, 
                                              levels = rev(c("Albacore", "Bluefin tuna", "Broadbill Swordfish",
                                      "Blue Shark","Shortfin Mako Shark", "Common Thresher Shark", "Bigeye Thresher Shark",
                                      "Short-Beaked Common Dolphin", "Long-Beaked Common Dolphin", "Northern Right Whale Dolphin")))
  # Plot
  # Histograms
  p1 <- ggplot(subset(hms_predator_lengths, !is.na(Predator_Length)), aes(y = Predator_Species, x = Predator_Length, 
                                                                          height = stat(density),fill = Predator_Species)) +
    geom_density_ridges(stat="binline", bins = 75, scale = 0.99, draw_baseline = FALSE, alpha = 0.75)+
    scale_fill_manual("Predator Species", values = c(paletteSp)) + 
    coord_cartesian(xlim = c(min(hms_predator_lengths$Predator_Length), 
                             max(hms_predator_lengths$Predator_Length)), expand = 0, clip = "off")+
    theme_classic(base_size = 14) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none",
           axis.text.x = element_text(angle = 45, vjust = 0.75))
 
   # Density plot
  p2 <- ggplot(subset(hms_predator_lengths, !is.na(Predator_Length)), 
               aes(y = Predator_Species, x = Predator_Length, fill = Predator_Species)) +
    geom_density_ridges(color = "black", scale = 0.9,alpha =0.5,quantile_lines = T,
                        quantile_fun = function(x,...)median(x), show.legend = T)+
    scale_fill_manual("Predator Species", values = c(paletteSp)) +
    scale_y_discrete(expand = c(0,0))+
    coord_cartesian(xlim = c(min(hms_predator_lengths$Predator_Length), 
                             max(hms_predator_lengths$Predator_Length)), expand = 0, clip = "off")+
    theme_classic(base_size = 14) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none",
           axis.text.x = element_text(angle = 45, vjust = 0.5))

   # Version of selected plot without labels for assembly in Photoshop
  p3 <- ggplot(subset(hms_predator_lengths, !is.na(Predator_Length)), aes(y = Predator_Species, x = Predator_Length, 
                                                                          height = stat(density),fill = Predator_Species)) +
    geom_density_ridges(stat="binline", bins = 75, scale = 0.99, draw_baseline = FALSE, alpha = 0.75)+
    scale_fill_manual("Predator Species", values = c(paletteSp)) + 
    coord_cartesian(xlim = c(min(hms_predator_lengths$Predator_Length), 
                             max(hms_predator_lengths$Predator_Length)), expand = 0, clip = "off")+
    theme_classic(base_size = 14) +
    theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none",
           axis.text = element_blank(), axis.title = element_blank())
  p3
  if(savePlot == "yes") {
    ggsave("./plots/predatorLengths_histogram.tiff", plot = p1, width = 3000, height = 2500, 
           units = "px", dpi = 400, compression = "lzw")
    ggsave("./plots/predatorLengths_density.tiff", plot = p2, width = 3000, height = 2500, 
           units = "px", dpi = 400, compression = "lzw")
    ggsave("./plots/predatorLengths_histogram_NoLabels.tiff", plot = p3, width = 3000, height = 2500, 
           units = "px", dpi = 400, compression = "lzw")
  }
  if(saveDF == "yes") {
    saveRDS(hms_predator_lengths, "./plots/dataFrames/plotData_plotLengthsPredators.rds")
  }
  return(p1) 
}
