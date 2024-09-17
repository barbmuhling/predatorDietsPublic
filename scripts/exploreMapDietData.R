#################################################################################################################
# Code to build Fig 1, showing predator locations and sizes
#################################################################################################################

library(ggplot2)
library(viridis)
library(scales)
library(tidyr)
library(dplyr)

# Map prep
coast <- borders(database = "world", fill = "grey90", colour = "black", xlim = c(-140, -100), ylim = c(0, 60))
# Palette
paletteSpecies <- data.frame("grp" = c("Albacore","Bluefin tuna", 
                                       "Long-Beaked Common Dolphin","Short-Beaked Common Dolphin", "Northern Right Whale Dolphin",
                                       "Common Thresher Shark","Blue Shark", "Shortfin Mako Shark", "Bigeye Thresher Shark",
                                       "Broadbill Swordfish"), 
                             "col" = c("#5073AB", "#090963","black","grey30", "grey70","#F2C2D2", "#CF977C", "#AB4467", "#8A401C",
                                       "#ACAFD2"))

# Load the datasets required: First prey composition
hms_prey <- read.csv("./data/hms_prey_withTraits.csv") 
# Then the file describing prey taxonomy
prey_taxonomy <- read.csv("./data/hms_prey_taxonomy.csv")
# Aggregate to number of gut samples per predator
predators <- data.frame(hms_prey %>% group_by(Predator_ID, Predator_Species, Year, Region) %>%
  summarise(Prey_N = sum(Prey_N)))
# Region as factor
predators$Region <- factor(predators$Region, levels = c("SoCal", "CenCal"))

#################################################################################################################
# Plots to show sampling coverage. Year/species with < 5 observations are shown in grey
spAgg <- aggregate(Predator_ID ~ Predator_Species + Year, predators, FUN = length)
# Re-order to group similar predators
spAgg$Predator_Species <- factor(spAgg$Predator_Species, levels = rev(c("Albacore", "Bluefin tuna", "Broadbill Swordfish", 
                                                                        "Blue Shark", "Shortfin Mako Shark", 
                                                                        "Common Thresher Shark", "Bigeye Thresher Shark",
                                    "Short-Beaked Common Dolphin", "Long-Beaked Common Dolphin", "Northern Right Whale Dolphin")))
# Plot
pSpYr <- ggplot() + geom_tile(data = subset(spAgg, Predator_ID > 2),
                              aes(x = Year, y = Predator_Species, fill = Predator_ID)) +
  scale_fill_gradientn(colours = c("#23045E", "#5092D9", "#C356DB", "#FF8FC3"),"No. \nSamples", limits = c(2, 100), oob = squish)+
  geom_tile(data = subset(spAgg, Predator_ID <= 2), aes(x = Year, y = Predator_Species), fill = "gray80") + 
  geom_vline(xintercept = 2015.5, linetype = "dashed")+
  theme_bw(base_size = 14) +
  scale_x_continuous(expand = c(0,0))+
  theme(legend.position = "bottom", panel.grid = element_blank())
# pSpYr
# Save
ggsave("./plots/sampleCoverageBySpeciesYear_v2018.tiff", plot = pSpYr, width = 3000, height = 2000, 
       units = "px", dpi = 400, compression = "lzw")

#################################################################################################################
# Maps of sampling coverage as hulls to preserve confidentiality
spp <- levels(spAgg$Predator_Species)
# Load processed hulls
allHulls <- readRDS("./data/sppHulls.rds")

# Sort out palette
paletteSpecies$grp <- factor(paletteSpecies$grp, levels = spp)
paletteSpecies <- paletteSpecies[order(paletteSpecies$grp),]
paletteSp <- paletteSpecies$col

# Map
pMap2 <- ggplot(data = allHulls) + 
  geom_path(aes(x = Longitude, y = Latitude, group = Predator_Species, color = Predator_Species), 
            linewidth = 1.75, na.rm = TRUE) +
  geom_polygon(aes(x = Longitude, y = Latitude, group = Predator_Species, fill = Predator_Species), color = NA,
  alpha = 0.5, na.rm = TRUE) +
  scale_color_manual("", values = rev(c(paletteSp))) + # guides(color = guide_legend(nrow = 4, byrow = TRUE)) +
  scale_fill_manual(values = rev(c(paletteSp))) + guides(fill = "none", color = "none") +
  coast + coord_quickmap(xlim = c(-130, -115), ylim = c(29.25, 41)) + scale_x_continuous(breaks = c(-128, -118)) + theme_bw() +
  facet_wrap(~ Predator_Species, ncol = 5) +
  theme(legend.position = "none",
        strip.background = element_blank(), strip.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
# pMap2
# Save
ggsave("./plots/sampleCoverageByRegionMap_polygons.tiff", plot = pMap2, 
       width = 4000, height = 1500, units = "px", dpi = 400, compression = "lzw")