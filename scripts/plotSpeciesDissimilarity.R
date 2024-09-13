#################################################################################################################
# Broad goal is to quantify how stable separation distance is relative to observed niche of each predator.
# See May and MacArthur a starting point https://www.pnas.org/doi/epdf/10.1073/pnas.69.5.1109.
# Standard deviation in dissimilarity between predators will be used as a proxy for separation distance. 
# Within group standard deviation dissimilarity will be used as a proxy for variability in observed niche for each predator. 
#################################################################################################################

plotSpeciesDissimilarity <- function(nmdsTS, savePlot) {
  toPerm <- nmdsTS$accumSub
  # Get species matrix
  spPerm <- toPerm[!names(toPerm) %in% c("gutsPerSpPerYr")]
  # Generate distance matrix
  dist <- vegdist(spPerm[!names(spPerm) %in% c("Predator_Species", "Year")], method = "horn")
  # Convert distance object to matrix
  distmat <- as.matrix(dist)
  # Add appropriate labels to rows/columns
  lab <- paste0(spPerm$Predator_Species, "_", spPerm$Year)
  dimnames(distmat) <- list(lab1 = lab, lab2 = lab) 
  
  # Melt longer, add labels
  distLong <- reshape2::melt(distmat, value.name = "distance")
  labs1 <- data.frame(str_split_fixed(distLong$lab1, "_", n = 2))
  colnames(labs1) <- c("species1", "year1")
  labs2 <- data.frame(str_split_fixed(distLong$lab2, "_", n = 2))
  colnames(labs2) <- c("species2", "year2")
  distLong <- cbind(distLong, labs1, labs2)
  # Years as numeric
  distLong$year1 <- as.numeric(distLong$year1)
  distLong$year2 <- as.numeric(distLong$year2)
  
  # Get all unique combinations of species
  spp <- unique(distLong$species1)
  sppAll <- data.frame(expand.grid("species1" = spp, "species2" = spp))
  # Remove species-species comparisons
  sppAll$remove <- ifelse(sppAll$species1 == sppAll$species2, "remove", "ok")
  sppAll <- subset(sppAll, remove == "ok")
  sppAll$remove <- NULL
 
  # Loop through this df to calculate within and between stdevs for all possible combinations of species
  # Species-species comparisons only within common years
  sdOut <- data.frame(matrix(nrow = nrow(sppAll), ncol = 9))
  colnames(sdOut) <- c("species1", "species2", "sdWithinSp1", "sdWithinSp2", "meanWithinSp1", "meanWithinSp2",
                       "sdAcrossBothSp", "meanAcrossBothSp","noOverlappingYrs")
  for(i in 1:nrow(sdOut)) {
    # Subset to just 2 species (e.g. albacore and blue shark)
    distSub <- subset(distLong, species1 == sppAll$species1[i] & species2 == sppAll$species2[i])
    # Add to df
    sdOut$species1[i] <- as.character(sppAll$species1[i])
    sdOut$species2[i] <- as.character(sppAll$species2[i])
    # Just for common years
    distSubYrs <- subset(distSub, year1 == year2) # e.g. albacore and blue shark have 5 years where both were sampled
    sdOut$noOverlappingYrs[i] <- nrow(distSubYrs)
    # If there are no years in common, skip 
    if(nrow(distSubYrs) == 0){
      sdOut[i,] <- NA 
      next
    }  
    # St dev of distance for these common years
    sdOut$meanAcrossBothSp[i] <- mean(distSubYrs$distance)
    # st dev of distance for these common years
    sdOut$sdAcrossBothSp[i] <- sd(distSubYrs$distance)
    # Now calculate within-species st devs, but only for these 5 shared years
    distSubSp1 <- subset(distLong, species1 == sppAll$species1[i] & species2 == sppAll$species1[i])
    # Subset to the relevant years
    distSubSp1 <- subset(distSubSp1, year1 %in% distSubYrs$year1 & year2 %in% distSubYrs$year1)
    # Remove same-same comparisons
    distSubSp1 <- subset(distSubSp1, distance > 0)
    # Calculate mean
    sdOut$meanWithinSp1[i] <- mean(distSubSp1$distance)
    # Calculate st dev
    sdOut$sdWithinSp1[i] <- sd(distSubSp1$distance)
    # Repeat for species 2
    distSubSp2 <- subset(distLong, species1 == sppAll$species2[i] & species2 == sppAll$species2[i])
    # Subset to the relevant years
    distSubSp2 <- subset(distSubSp2, year1 %in% distSubYrs$year1 & year2 %in% distSubYrs$year1)
    # Remove same-same comparisons
    distSubSp2 <- subset(distSubSp2, distance > 0)
    # Calculate mean
    sdOut$meanWithinSp2[i] <- mean(distSubSp2$distance)
    # Calculate st dev
    sdOut$sdWithinSp2[i] <- sd(distSubSp2$distance)
  }
  
  colnames(sdOut) <- c("species1", "species2", "sdWithinSp1", "sdWithinSp2", "meanWithinSp1", "meanWithinSp2", 
                       "sdAcrossBothSp","meanAcrossBothSp", "noOverlappingYrs")
  
  # Reshape so that we're isolating each species (whether it was sp1 or sp 2)
  sdOutLong <- pivot_longer(sdOut[c("species1", "species2", "sdWithinSp1",  "sdAcrossBothSp")], 
                    cols = c( "sdWithinSp1","sdAcrossBothSp"), 
                    names_to = "sdType", values_to = "sd")
  # Reshape so that we're  keeping shape for sign test
  sdOutCrop <- sdOut[,c("species1", "species2", "sdWithinSp1", "sdAcrossBothSp")]
  colnames(sdOutCrop) <- c("species1", "species2", "sdWithin", "sdAcross")
  sdOutLong$spType <- paste0(sdOutLong$species1, "_", sdOutLong$sdType) 
  sdOutLong$species2 <- ifelse(sdOutLong$sdType == "sdWithinSp1",sdOutLong$species1,  sdOutLong$species2)
  # Define which predators are present here
  paletteSpecies <- subset(paletteSpecies, grp %in% sdOutLong$species1) 
  # Order predators by schoeners clusters
  schoeners_order <- c("Common Thresher Shark","Long-Beaked Common Dolphin",
                       "Bluefin tuna", "Albacore",
                       "Blue Shark", "Broadbill Swordfish","Shortfin Mako Shark",
                       "Short-Beaked Common Dolphin")
  paletteSpecies$grp <- reorder.factor(paletteSpecies$grp, new.order = schoeners_order)
  paletteSpecies <-paletteSpecies[order(paletteSpecies$grp),]
  paletteSp <- paletteSpecies$col # Groups predators by color 
  sdOutLong$species1 <- reorder.factor(sdOutLong$species1, new.order = schoeners_order)
  sdOutLong$species2 <- reorder.factor(sdOutLong$species2, new.order = schoeners_order)
  sdOut$species1 <- reorder.factor(sdOut$species1, new.order = schoeners_order)
  sdOut$species2 <- reorder.factor(sdOut$species2, new.order = schoeners_order)
  
  paletteSpMOD <- c(paletteSp[4], "white", paletteSp[5], "white", paletteSp[3], "white", paletteSp[6], "white", 
                    paletteSp[1], "white", paletteSp[2], "white", paletteSp[8], "white", paletteSp[7], "white")
  
  p1 <- ggplot(sdOutLong, aes(x = species1, y = sd, group = spType, fill = spType, color = species1)) + 
    geom_boxplot(alpha = 0.5, size =1.05)+theme_bw()+
    scale_fill_manual(values = paletteSpMOD, guide = NULL)+
    scale_color_manual(values = paletteSp, guide = NULL)+
    scale_y_continuous(limits = c(0.05, 0.35))+
    theme(panel.grid.major.x = element_blank())
  p1
  
  p2 <- ggplot(subset(sdOutLong, sdType == "sdAcrossBothSp"), aes(x = species1, y = sd)) +
    geom_point(aes(color = species2),size = 2.5, alpha = 0.5)+
    geom_point(aes(color = species2),shape = 1,size = 2.5, stroke =1.5)+theme_bw()+
    scale_color_manual(values = paletteSp, guide = NULL)+
    scale_fill_manual(values = paletteSp, guide = NULL)+
    scale_y_continuous(limits = c(0.05, 0.35))+
    theme(panel.grid = element_blank())
  p2
  
  p3 <- ggplot(sdOut, aes(x = meanWithinSp1, y = sdWithinSp1)) +
    geom_point(aes(color = species1),size = 4, alpha = 0.75) +
    geom_smooth(method = "lm", se=T, color = "black") +
    scale_color_manual(values = paletteSp, guide = NULL)+
    facet_wrap(~ species1)+
    labs(x = "dissimilarity", y = "standard deviation of dissimilarity")+
    theme(panel.grid = element_blank()) + theme_bw()
  p3
  
  p3 <- p3 +
    geom_point(aes(x = meanAcrossBothSp, y = sdAcrossBothSp,color = species2),shape = 15, size = 4, alpha = 0.75)+
    geom_smooth(method = "lm",se=T, aes(x = meanAcrossBothSp, y = sdAcrossBothSp), color = "grey")+
    scale_color_manual(values = paletteSp, guide = NULL)+
    facet_wrap(~ species1)+
    theme(panel.grid = element_blank())+theme_bw()
  p3
  
  # Loop through this df to calculate within and between stdevs for all possible combinations of species
  # Species-species comparisons only within common years
  WicoxonRankOut <- data.frame(matrix(nrow = length(unique(sdOutCrop$species1)), ncol = 3))
  # Running wilcoxon rank sum test test for a greater than b AND a less than b hypotheses.
  colnames(WicoxonRankOut) <- c("species1","p_value_greater", "p_value_less")
  for(i in 1:nrow(WicoxonRankOut)) {
    testSub <- subset(sdOutCrop, species1 == sdOutCrop$species1[i])
    # Add to df
    WicoxonRankOut$species1[i] <- as.character(testSub$species1[1])
    wtestgreat <- wilcox.test(x=testSub$sdWithin,y =testSub$sdAcross, paired = T, alternative = "greater")
    WicoxonRankOut$p_value_greater[i] <- wtestgreat$p.value
    wtestless <- wilcox.test(x=testSub$sdWithin,y =testSub$sdAcross, paired = T,  alternative = "less")
    WicoxonRankOut$p_value_less[i] <- wtestless$p.value
  }

  # Save
  if(savePlot == "yes") {
    ggsave("./plots/ordinations/Boxplots_withinVSacross_dissimilarity_SD.tiff",  plot = p1, width = 4000, 
           height = 2000, units = "px", dpi = 400, compression = "lzw")
    ggsave("./plots/ordinations/Points_across_dissimilarity_SD.tiff",  plot = p2, width = 2500, 
           height = 2000, units = "px", dpi = 400, compression = "lzw")
  }
  # Return a list of plots and outputs
  out <- list("speciesBoxPlot" = p1,"AcrossPointPlot" = p1, "dissims" = sdOutLong, "wilcoxonOut"=WicoxonRankOut)
  return(out) 
}