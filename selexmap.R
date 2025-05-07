# install package for interfacing with Cytoscape
if(!require("BiocManager")) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install("RCy3", force=TRUE)
}

# install mappeR package
if(!require("mappeR")) {
  install.packages("mappeR")
}

# load required packages
library(mappeR)
library(readxl)
library(stringdist)
library(usedist)

source("cytoscape-stuff.R")

# import from excel file
selex_data = na.omit(as.data.frame(read_xlsx("data/selex_data.xlsx", sheet = "Data")))

# individual observations are aptamers
rownames(selex_data) = selex_data$Name

# define a function to get the edit distance between two aptamers
get_edit_distance <- function(apt1, apt2) {
  return(stringdist(apt1["Aptamer"], apt2["Aptamer"], method = "lv"))
}

# generate distance matrix for aptamers using edit distance
edit_dists = dist_make(as.matrix(selex_data), get_edit_distance)

# create ball mapper object with edit distances and specified epsilon value
eps = 2
aptamerballs = create_ball_mapper_object(selex_data, edit_dists, eps)

# get aptamers per vertex of the ballmapper graph
balled_aptamers = lapply(aptamerballs[[1]]$data, function(x) selex_data[unlist(strsplit(x, ", ")),])

# calculate median round amounts for the aptamers in each vertex
rounds = c("RP10M 0A", "RP10M 0B", "RP10M 0C", "RP10M 1A", "RP10M 1B", "RP10M 1C", "RP10M 2A", "RP10M 2B", "RP10M 2C", "RP10M 3A", "RP10M 3B", "RP10M 3C", "RP10M 4", "RP10M 5", "RP10M 6", "RP10M 7", "RP10M 8", "RP10M 9")
round_averages = sapply(balled_aptamers, function(x) sapply(rounds, function(r) mean(selex_data[x$Name, r])))

# calculate median human affinity values aptamers per vertex
median_human_affinity = sapply(balled_aptamers, function(x) median(selex_data[x$Name, "ASSET hVSMC:hEC"]))
mean_human_affinity = sapply(balled_aptamers, function(x) mean(selex_data[x$Name, "ASSET hVSMC:hEC"]))
max_human_affinity = sapply(balled_aptamers, function(x) max(selex_data[x$Name, "ASSET hVSMC:hEC"]))

# calculate median mouse affinity values of aptamers per vertex
mean_human_affinity = sapply(balled_aptamers, function(x) mean(selex_data[x$Name, "ASSET mVSMC:mEC"]))
mean_mouse_affinity = sapply(balled_aptamers, function(x) mean(selex_data[x$Name, "ASSET mVSMC:mEC"]))
max_mouse_affinity = sapply(balled_aptamers, function(x) max(selex_data[x$Name, "ASSET mVSMC:mEC"]))

# add affinity info to mapper graph data
aptamerballs[[1]][rounds] = t(round_averages)
aptamerballs[[1]]$mean_human_affinity = log2(mean_human_affinity)
aptamerballs[[1]]$mean_mouse_affinity = log2(mean_mouse_affinity)
aptamerballs[[1]]$enrichment = log2(aptamerballs[[1]]$"RP10M 9"/aptamerballs[[1]]$"RP10M 3A")

# send mapper data to cytoscape
cymapper(aptamerballs, is_ballmapper = TRUE)

assets = selex_data$`ASSET hVSMC:hEC`
num_patches = 5
percent_overlap = 20

assetmapper = create_1D_mapper_object(selex_data, edit_dists, selex_data$`ASSET hVSMC:hEC`, create_width_balanced_cover(min(assets), max(assets), num_patches, percent_overlap))

cymapper(assetmapper)

