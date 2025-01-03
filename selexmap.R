library(xlsx)
library(usedist)
library(stringdist)
library(devtools)

install.packages("BiocManager", repos = "https://cloud.r-project.org")

BiocManager::install("RCy3")

if(!require("mappeR")) {
  install_github("https://github.com/Uiowa-Applied-Topology/mappeR/tree/cytoscape")
}

library(mappeR)

selex_data = na.omit(as.data.frame(read.xlsx("data/selex_data.xlsx", sheetName = "Data")))

rownames(selex_data) = selex_data$Name

get_edit_distance <- function(apt1, apt2) {
  return(stringdist(apt1["Variable"], apt2["Variable"], method = "lv"))
}

edit_dists = dist_make(as.matrix(selex_data), get_edit_distance)

aptamerballs = create_ball_mapper_object(selex_data, edit_dists, 6)

balled_aptamers = lapply(aptamerballs[[1]]$data, function(x) selex_data[unlist(strsplit(x, ", ")),])

mean_human_affinity = sapply(balled_aptamers, function(x) mean(selex_data[x$Name, "ASSET.hVSMC.hEC"]))
mean_mouse_affinity = sapply(balled_aptamers, function(x) mean(selex_data[x$Name, "ASSET.mVSMC.mEC"]))

aptamerballs[[1]]$mean_human_affinity = mean_human_affinity
aptamerballs[[1]]$mean_mouse_affinity = mean_mouse_affinity

cymapper(aptamerballs, is_ballmapper = TRUE)
