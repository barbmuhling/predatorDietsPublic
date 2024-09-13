###############################################################################################################
# Main script to call functions to process data, complete multivariate analyses,
# plot ordinations and examine foraging depths.
# Originally written by Elan Portner and Barb Muhling, maintained by Barb Muhling 
# (barbara.muhling@noaa.gov)
###############################################################################################################

###############################################################################################################
########################################## Load package libraries #############################################
###############################################################################################################
library(labdsv);library(vegan); library(plyr); library(dplyr);library(ggplot2);
library(tidyr); library(ggrepel); library(viridis); library(scales); library(readxl);library(stringr);
library(dendextend); library(lubridate);
library(ggridges); library(gdata); library(ggpubr); library(iNEXT); library(grid)
# install pairwise adonis package from github if needed
# library(remotes); library(devtools)
# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

################################################################################################################
####################################### Load required data files ###############################################
################################################################################################################
hms_prey <- read.csv("./data/hms_prey_withTraits.csv") 
hms_prey_lengths <- read.csv("./data/hms_prey_lengths_subset.csv")
hms_predator_lengths <- read.csv("./data/hms_predator_lengths_subset.csv")
hms_predator_depths <- read.csv("./data/hms_predator_depths.csv")
prey_taxonomy <- read.csv("./data/hms_prey_taxonomy.csv")

##############################################################################################################
########################################## Define color palettes #############################################
##############################################################################################################
# Predator species palette
paletteSpecies <- data.frame("grp" = c("Albacore","Bluefin tuna", 
                                     "Long-Beaked Common Dolphin","Short-Beaked Common Dolphin", 
                                     "Northern Right Whale Dolphin", "Common Thresher Shark","Blue Shark", 
                                     "Shortfin Mako Shark", "Bigeye Thresher Shark", "Broadbill Swordfish"), 
                             "col" = c("#5073AB", "#090963", "black", "grey30", "grey70", "#F2C2D2", "#CF977C", 
                                       "#AB4467", "#8A401C", "#ACAFD2"))
# Prey species palette
palettePrey <- data.frame("Prey_LPT" = c("Bathylagidae","Carangidae","Clupeidae","Engraulidae","Merlucciidae",
                                     "Microstomatidae","Mugilidae" ,"Myctophidae","Paralepididae","Scomberesocidae",
                                     "Scombridae","Scopelarchidae","Sebastidae", "Synodontidae","Trachipteridae",   
                                     "Alloposidae","Argonautidae","Enoploteuthidae","Gonatidae","Histioteuthidae",
                                     "Loliginidae","Octopodidae" ,"Octopoteuthidae","Ommastrephidae","Onychoteuthidae",
                                     "Munididae"),
                             "col" = c("#CEE6F2","#88B7CF","#508DAB","#266787","#094563",
                                    "#CEDCF2","#88A3CF","#5073AB","#264B87","#092C63",
                                    "#CECEF2","#8888CF","#5050AB","#262687","#090963",
                                    "#F2C2D2","#CF7C97","#AB4467","#8A1C40","#660022",
                                    "#E6BFAC","#CF977C","#AB6744","#8A401C","#662200","grey90"))

##############################################################################################################
################################### Compare prey lengths across predator taxa ################################
##############################################################################################################
source("./scripts/plotLengths.R")
pLengths <- plotLengths(hms_prey_lengths, paletteSpecies = paletteSpecies, savePlot = "yes", saveDF = "no")
# Call desired plot from list
pLengths$MaxLengthsThrough2015
# Stats of median prey size changes for PBF and SWO 1998-2015 vs 2016-2018. Measured items only
aggregate(Prey_Length ~ subset, pLengths$BluefinTemporal$data, FUN = median, na.rm = TRUE) 
aggregate(Prey_Length ~ subset, pLengths$SwordfishTemporal$data, FUN = median, na.rm = TRUE) 
# Mann-Whitney u-tests
wilcox.test(Prey_Length ~ subset, data = pLengths$BluefinTemporal$data) 
wilcox.test(Prey_Length ~ subset, data = pLengths$SwordfishTemporal$data) 

###############################################################################################################
####################################### Compare predator lengths ##############################################
###############################################################################################################
source("./scripts/plotLengthsPredators.R")
pLengths2 <- plotLengthsPredators(hms_predator_lengths, paletteSpecies = paletteSpecies, savePlot = "yes", saveDF = "no")
pLengths2

###############################################################################################################
################################# Set parameters for multivariate analyses ####################################
###############################################################################################################
# For creating the species matrix for NMDS
startYr <- 1998
endYr <- 2018
preyCutoff <- 0.025 # 0.025 means include prey species that contribute at least 2.5% in at least 1 predator
saveFile <- "yes" # Save the matrix to disk as a RDS
familyAgg <- "yes" # Aggregate to family level before building matrix for NMDS 
yrAgg <- "no" # average (aggregate) data across years? "yes" or "no"
k <- 4
try <- 20
trymax <- 100
saveNMDS <- "yes" # "yes" or "no"

################################################################################################################
########################################### Create species matrix ##############################################
################################################################################################################
# Create a matrix for multivariate analyses
source("./scripts/getSpMatrix.R")
# 1998 - 2015
accum <- getSpMatrix(hms_prey = hms_prey, prey_taxonomy = prey_taxonomy, startYr = startYr, endYr = endYr, endYrAc = 2015,
                     preyCutoff = preyCutoff, saveFile = saveFile, familyAgg = familyAgg)

# 1998 - 2018
accumfull <- getSpMatrix(hms_prey = hms_prey, prey_taxonomy = prey_taxonomy, startYr = startYr, endYr = endYr, endYrAc = endYr,
                         preyCutoff = preyCutoff, saveFile = saveFile, familyAgg = familyAgg)

################################################################################################################
################################## Plot Diversity Accummulation Curves #########################################
################################################################################################################
source("./scripts/plotAccumulationCurves.R")
# Can run curves by year per species, or just by species
type = "total" # can be "total" (all predator species on one plot) or "bySpecies" (plots by year per predator) or "both" 
# Note that short-beaked common dolphin accum curves bySpecies "break" when nboot is too high. 
pCurves <- plotAccumulationCurves(accum = accumfull, paletteSpecies = paletteSpecies, type = type, 
                                  savePlot = "yes", nboot = 150, q = 0) 

################################################################################################################
################################################# Run NMDS #####################################################
################################################################################################################
# Using specified subset of the full "accum" matrix
source("./scripts/runNMDS.R")
# Run function to output a list with the NMDS object, and the data matrix used to create it
nmds <- runNMDS(accum = accum, yrAgg = yrAgg, k = k, try = try, trymax = trymax, saveNMDS = "yes", 
                preyCutoff = preyCutoff, familyAgg = familyAgg, regionAgg = "no")
# # Can also load an NMDS saved previously
# nmds <- readRDS(file = paste0("./data/nmds_", preyCutoff, "_", "familyAgg", familyAgg, "_yrAgg", yrAgg,
#                       "_", startYr, "_", endYr, ".rds")) 

################################################################################################################
############################################# Plot NMDS results ################################################
################################################################################################################
# Generate plots visualizing NMDS results. Plots are saved if specified
source("./scripts/plotScatter.R")
# Scatterplot
pScatter <- plotScatter(myNMDS = nmds$nmds, myAccum = nmds$accumSub, preyCutoff = preyCutoff,  
                        familyAgg = familyAgg, yrAgg = yrAgg, paletteSpecies = paletteSpecies, 
                        savePlot = "yes", saveDF = "no")
# Returns a list of several plots: can visualize whichever one you like
pScatter$facetNMDS

source("./scripts/plotScatterVectors.R")
# Scatterplot with vectors. Note hyp argument to select taxa for which to plot vectors
pScatterVectors <- plotScatterVectors(myNMDS = nmds$nmds, myAccum = nmds$accumSub, preyCutoff = preyCutoff, 
                                     familyAgg = familyAgg, yrAgg = yrAgg, 
                                     paletteSpecies = paletteSpecies, savePlot = "yes", 
                                     saveDF = "no", hypCut = 0.25)
pScatterVectors$vectorLabels
pScatterVectors$vectors

###############################################################################################################
############################################### Time series analyses ##########################################
###############################################################################################################
source("./scripts/plotOrdinationTrajectory.R")
# All predators except northern right whale dolphin and bigeye thresher shark
# returns a list of plots and summary stats
pTrajectory <- plotOrdinationTrajectory(hms_prey = hms_prey, prey_taxonomy = prey_taxonomy, startYrTS = startYr, endYrTS = endYr, 
                                         preyCutoff = preyCutoff, familyAggTS = "yes", paletteSpecies = paletteSpecies,
                                         spp = c("Albacore", "Blue Shark", "Bluefin tuna", 
                                                 "Broadbill Swordfish", "Common Thresher Shark", "Long-Beaked Common Dolphin", 
                                                 "Short-Beaked Common Dolphin", "Shortfin Mako Shark"),
                                         savePlot = "yes", saveDF = "no")
pTrajectory$TrajectoriesByYear # Or whichever plot you want to visualize

# Calculate and plot mean dissimilarity between species (common years), and within species
source("./scripts/plotSpeciesDissimilarity.R")
pDissim <- plotSpeciesDissimilarity(nmdsTS = pTrajectory$nmdsTS, savePlot = "yes") # Returns dfs and plots
pDissim$speciesBoxPlot

###############################################################################################################
######################## Schoener's Index of overlap, using prey matrix calculated above ######################
###############################################################################################################
source("./scripts/calculateSchoeners.R")
pSchoen <- calculateSchoeners(hms_prey = hms_prey, accum = accum, savePlot = "yes", saveDF = "no")
pSchoen$schCluster

###############################################################################################################
########################################### Multivariate Statistics ###########################################
###############################################################################################################
source("./scripts/runMultivariateStats.R") # slow, lower nPerm allows testing
runMultivariateStats(accumSub = nmds$accumSub, nPerm = 999, saveOutputs = "no")

###############################################################################################################
##################################### Prey family barplots by predator ########################################
###############################################################################################################
source("./scripts/plotFamilyLevelDietBars.R")
familyBarPlot <- plotDietBars(preyCutoff = preyCutoff, startYr = startYr, endYr = endYr, hms_prey = hms_prey, 
                              cropStart = startYr, cropEnd = 2015,
                                     prey_taxonomy = prey_taxonomy, savePlot = "yes", saveDF = "no",
                              spp = c("Albacore","Bigeye Thresher Shark" ,"Blue Shark", "Bluefin tuna", 
                                      "Broadbill Swordfish", "Common Thresher Shark", "Long-Beaked Common Dolphin", 
                                      "Short-Beaked Common Dolphin", "Shortfin Mako Shark", "Northern Right Whale Dolphin"))
familyBarPlot$MeanDietNoLabs

# Generate temporal subsets for bluefin and swordfish
# Earlier time-period
familyBarPlotCROP <- plotDietBars(preyCutoff = preyCutoff, startYr = startYr, endYr = endYr, hms_prey = hms_prey, 
                              cropStart = 2006, cropEnd = 2015,
                              prey_taxonomy = prey_taxonomy, savePlot = "yes", saveDF = "no",
                              spp = c("Bluefin tuna", "Broadbill Swordfish"))
familyBarPlotCROP$MeanDietNoLabs

# Later time-period
familyBarPlotCROP2 <- plotDietBars(preyCutoff = preyCutoff, startYr = startYr, endYr = endYr, hms_prey = hms_prey, 
                                  cropStart = 2016, cropEnd = 2018,
                                  prey_taxonomy = prey_taxonomy, savePlot = "yes", saveDF = "no",
                                  spp = c("Bluefin tuna","Broadbill Swordfish"))
familyBarPlotCROP2$MeanDietNoLabs

#################################################################################################################
############################### Prey vertical habitat barplots by species/year ##################################
#################################################################################################################
source("./scripts/plotVerticalHabitats.R")
traitBarPlot <- plotVerticalHabitats(preyCutoff = 0, startYr = startYr, endYr = endYr, hms_prey = hms_prey, 
                                     cropStart = startYr, cropEnd = 2015,
                                     prey_taxonomy = prey_taxonomy, savePlot = "yes", saveDF = "no",
                                     spp = c("Albacore","Bigeye Thresher Shark" ,"Blue Shark", "Bluefin tuna", 
                                       "Broadbill Swordfish", "Common Thresher Shark", "Long-Beaked Common Dolphin", 
                                       "Short-Beaked Common Dolphin", "Shortfin Mako Shark", "Northern Right Whale Dolphin"))
traitBarPlot$horizontal
# Plot by time-period for bluefin and swordfish
traitBarPlotCROP <- plotVerticalHabitats(preyCutoff = 0, startYr = startYr, endYr = endYr, hms_prey = hms_prey, 
                                                cropStart = 2006, cropEnd = 2015,
                                                prey_taxonomy = prey_taxonomy, savePlot = "yes", saveDF = "no",
                                                spp = c("Bluefin tuna","Broadbill Swordfish"))
traitBarPlotCROP$horizontal

#################################################################################################################
############################ Prey inshore/offshore habitat barplots by species/year #############################
#################################################################################################################
source("./scripts/plotInshoreOffshoreHabitats.R")
inOffBarPlot <- plotInshoreOffshoreHabitats(preyCutoff = 0, startYr = startYr, endYr = endYr, hms_prey = hms_prey,
                                            cropStart = startYr, cropEnd = 2015,
                                     prey_taxonomy = prey_taxonomy, savePlot = "yes", saveDF = "no",
                                     spp = c("Albacore","Bigeye Thresher Shark" ,"Blue Shark", "Bluefin tuna", 
                                          "Broadbill Swordfish", "Common Thresher Shark", "Long-Beaked Common Dolphin", 
                                          "Short-Beaked Common Dolphin", "Shortfin Mako Shark", "Northern Right Whale Dolphin"))
inOffBarPlot$horizontal
# Plot by time-period for bluefin and swordfish
inOffBarPlotCROP <- plotInshoreOffshoreHabitats(preyCutoff = 0, startYr = startYr, endYr = endYr, hms_prey = hms_prey, 
                                         cropStart = 2006, cropEnd = 2015,
                                         prey_taxonomy = prey_taxonomy, savePlot = "yes", saveDF = "no",
                                         spp = c("Bluefin tuna","Broadbill Swordfish"))
inOffBarPlotCROP$horizontal

################################################################################################################
####################################### Calculate and plot depth of foraging ###################################
################################################################################################################
# Generate species matrix not aggregated to family
source("./scripts/getSpMatrix.R")
accumDepth <- getSpMatrix(hms_prey = hms_prey, prey_taxonomy = prey_taxonomy, startYr = startYr, endYr = endYr,
                          endYrAc = 2015, preyCutoff = preyCutoff, saveFile = "no", familyAgg = "no")
# Plot depth of foraging
source("./scripts/plotDepthForaging.R") 
# Note preyMetric can be "occurrence" or "abundance"
depthForagingPlot <- plotDepthForaging(accumDepth = accumDepth, hms_prey = hms_prey, preyMetric = "abundance",
                                       hms_predator_depths = hms_predator_depths, paletteSpecies = paletteSpecies, 
                                       savePlot = "yes", saveDF = "no") 
depthForagingPlot
