if(!require("BiocManager")) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install("RCy3", force=TRUE)
}

if(!require("devtools")) {
  install.packages("devtools")
}

if(!require("mappeR")) {
  install_github("https://github.com/Uiowa-Applied-Topology/mappeR/tree/cytoscape")
}

library(devtools)
library(mappeR)
library(readxl)
library(stringdist)
library(usedist)

# import from excel file
selex_data = na.omit(as.data.frame(read_xlsx("data/selex_data.xlsx", sheet = "Data")))

# observations are aptamers
rownames(selex_data) = selex_data$Name

# get the edit distance between two aptamers
get_edit_distance <- function(apt1, apt2) {
  return(stringdist(apt1["Variable"], apt2["Variable"], method = "lv"))
}

# generate distance matrix for aptamers using edit distance
edit_dists = dist_make(as.matrix(selex_data), get_edit_distance)

# create ball mapper object with edit distances
aptamerballs = create_ball_mapper_object(selex_data, edit_dists, 10)

# get aptamers in vertices of ballmapper graph
balled_aptamers = lapply(aptamerballs[[1]]$data, function(x) selex_data[unlist(strsplit(x, ", ")),])

# find mean human/mouse affinity per ballmapper vertex
mean_human_affinity = sapply(balled_aptamers, function(x) mean(selex_data[x$Name, "ASSET.hVSMC.hEC"]))
mean_mouse_affinity = sapply(balled_aptamers, function(x) mean(selex_data[x$Name, "ASSET.mVSMC.mEC"]))

# add mean affinity to vertex data
aptamerballs[[1]]$mean_human_affinity = mean_human_affinity
aptamerballs[[1]]$mean_mouse_affinity = mean_mouse_affinity

# send ballmapper data to cytoscape
cymapper(aptamerballs, is_ballmapper = TRUE)

