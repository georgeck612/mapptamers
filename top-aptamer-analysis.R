if(!require("BiocManager")) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install("RCy3", force=TRUE)
}

library(mappeR)
library(RCy3)


# data input and distance matrix creation ---------------------------------


# I don't think there are actually missing values but I can't trust myself
node_data = na.omit(as.data.frame(read.csv("data/top-nodes.csv")))
edge_data = na.omit(as.data.frame(read.csv("data/top-edges.csv")))

node_data[,c("hVSMC.hEC", "Log.10.RP10M.9", "Log2.hVSMC.hEC", "Log2.mVSMC.mEC", "Log2.R3.9", "mVSMC.mEC")] = round(node_data[,c("hVSMC.hEC", "Log.10.RP10M.9", "Log2.hVSMC.hEC", "Log2.mVSMC.mEC", "Log2.R3.9", "mVSMC.mEC")], 3)

# making sure R doesn't have any cause to yell at me
row.names(node_data) = node_data$name

total_aptamers = nrow(node_data)

# create empty edit distance matrix
edit_dists = matrix(0, nrow = total_aptamers, ncol = total_aptamers)

# flatten interactions column into a list of strings formatted as
# ["source aptamer", "to", "target aptamer"]
squished_edges = unlist(strsplit(edge_data$name, " "))

# recover sources/targets via above formatting
sources = squished_edges[seq(1, length(squished_edges), 3)]
targets = squished_edges[seq(3, length(squished_edges), 3)]

# grab entire source/target aptamer locations in mother dataset
row_indices = match(sources, node_data$name)
col_indices = match(targets, node_data$name)

# populate edit distance matrix
edit_dists[cbind(row_indices, col_indices)] = edge_data$editDistance
edit_dists[cbind(col_indices, row_indices)] = edge_data$editDistance

# not sure if necessary but I did it last time and it doesn't hurt to leave it
row.names(edit_dists) = node_data$name
colnames(edit_dists) = node_data$name

# create empty tree distance matrix
tree_dists = matrix(0, nrow = total_aptamers, ncol = total_aptamers)

# populate tree distance matrix
tree_dists[cbind(row_indices, col_indices)] = edge_data$treeDistance
tree_dists[cbind(col_indices, row_indices)] = edge_data$treeDistance

# again I'm being paranoid
row.names(tree_dists) = node_data$name
colnames(tree_dists) = node_data$name


# mapper graph creation ---------------------------------------------------

#
# # function which makes a ballmapper graph and populates it with data
# create_ballmapptamer_graph <- function(dists, eps) {
#   # mapper time
#   mapptamer = create_ball_mapper_object(node_data, dists, eps)
#
#   # get aptamers in vertices of mapper graph
#   aptamer_balls = lapply(mapptamer[[1]]$data, function(x) node_data[unlist(strsplit(x, ", ")),])
#
#   # calculate median selex data across mapper vertices
#   median_log10_final_selex_read = sapply(aptamer_balls, function(x) median(node_data[x$name, "Log.10.RP10M.9"]))
#   median_log2_selex_enrichment = sapply(aptamer_balls, function(x) median(node_data[x$name, "Log2.R3.9"]))
#
#   # calculate median ASSET data across mapper vertices
#   median_human_affinity = sapply(aptamer_balls, function(x) median(node_data[x$name, "hVSMC.hEC"]))
#   median_log2_human_affinity = sapply(aptamer_balls, function(x) median(node_data[x$name, "Log2.hVSMC.hEC"]))
#   median_mouse_affinity = sapply(aptamer_balls, function(x) median(node_data[x$name, "mVSMC.mEC"]))
#   median_log2_mouse_affinity = sapply(aptamer_balls, function(x) median(node_data[x$name, "mVSMC.mEC"]))
#
#   # attach calculated info to mapper dataframe
#   mapptamer[[1]]$median_log10_final_selex_read = median_log10_final_selex_read
#   mapptamer[[1]]$median_log2_selex_enrichment = median_log2_selex_enrichment
#   mapptamer[[1]]$median_human_affinity = median_human_affinity
#   mapptamer[[1]]$median_log2_human_affinity = median_log2_human_affinity
#   mapptamer[[1]]$median_mouse_affinity = median_mouse_affinity
#   mapptamer[[1]]$median_log2__mouse_affinity = median_log2_mouse_affinity
#
#   return(mapptamer)
# }

# function which makes a ballmapper graph and populates it with data
create_mapptamer_graph <- function(dists, filtered, cover, linkage_method) {
  # mapper time
  mapptamer = create_1D_mapper_object(node_data, dists, filtered, cover, clusterer = global_tallest_hierarchical_clusterer(linkage_method, as.dist(dists)))

  # get aptamers in vertices of mapper graph
  aptamer_balls = lapply(mapptamer[[1]]$data, function(x) node_data[unlist(strsplit(x, ", ")),])

  # calculate median selex data across mapper vertices
  median_log10_final_selex_read = sapply(aptamer_balls, function(x) median(node_data[x$name, "Log.10.RP10M.9"]))
  median_log2_selex_enrichment = sapply(aptamer_balls, function(x) median(node_data[x$name, "Log2.R3.9"]))

  # calculate median ASSET data across mapper vertices
  median_human_affinity = sapply(aptamer_balls, function(x) median(node_data[x$name, "hVSMC.hEC"]))
  median_log2_human_affinity = sapply(aptamer_balls, function(x) median(node_data[x$name, "Log2.hVSMC.hEC"]))
  median_mouse_affinity = sapply(aptamer_balls, function(x) median(node_data[x$name, "mVSMC.mEC"]))
  median_log2_mouse_affinity = sapply(aptamer_balls, function(x) median(node_data[x$name, "mVSMC.mEC"]))

  # attach calculated info to mapper dataframe
  mapptamer[[1]]$median_log10_final_selex_read = round(median_log10_final_selex_read, 3)
  mapptamer[[1]]$median_log2_selex_enrichment = round(median_log2_selex_enrichment, 3)
  mapptamer[[1]]$median_human_affinity = round(median_human_affinity, 3)
  mapptamer[[1]]$median_log2_human_affinity = round(median_log2_human_affinity, 3)
  mapptamer[[1]]$median_mouse_affinity = round(median_mouse_affinity, 3)
  mapptamer[[1]]$median_log2__mouse_affinity = round(median_log2_mouse_affinity, 3)

  return(mapptamer)
}

source("cytoscape-stuff.R")


# # this will crash when there are no edges, that's not a mappeR thing, that's an RCy3 thing
# for (eps in 2:8) {
#   mapptamer = create_ballmapptamer_graph(edit_dists, eps)
#   cymapper(mapptamer, mapptamer[[1]]$median_log10_final_selex_read, mapptamer[[1]]$median_human_affinity, mapptamer[[1]]$cluster_size)
# }
#
# for (eps in 2:39) {
#   mapptamer = create_ballmapptamer_graph(tree_dists, eps)
#   cymapper(mapptamer, mapptamer[[1]]$median_log10_final_selex_read, mapptamer[[1]]$median_human_affinity, mapptamer[[1]]$cluster_size)
# }

# NODE SIZE: selex enrichment
# NODE BORDER COLOR: human affinity
# NODE FILL COLOR: average distance to medoid
cytoviz1 = function(mapptamer, name) {
  size_col = "median_log2_selex_enrichment"
  size_data = sort(mapptamer[[1]][, size_col])
  border_color_col = "median_log2_human_affinity"
  border_color_data = sort(mapptamer[[1]][, border_color_col])
  fill_color_col = "tightness"
  fill_color_data = sort(mapptamer[[1]][, fill_color_col])
  visualize_mapper_data(mapptamer, size_col, size_data, border_color_col, border_color_data, fill_color_col, fill_color_data, name)
}

# NODE SIZE: human affinity
# NODE BORDER COLOR: selex enrichment
# NODE FILL COLOR: average distance to medoid
cytoviz2 = function(mapptamer, name) {
  size_col = "median_log2_human_affinity"
  size_data = sort(mapptamer[[1]][, size_col])
  border_color_col = "median_log2_selex_enrichment"
  border_color_data = sort(mapptamer[[1]][, border_color_col])
  fill_color_col = "tightness"
  fill_color_data = sort(mapptamer[[1]][, fill_color_col])
  visualize_mapper_data(mapptamer, size_col, size_data, border_color_col, border_color_data, fill_color_col, fill_color_data, name)
}

source("aptamer_clusterer.R")

percent_overlap = 25
affinity_projection = node_data$Log2.hVSMC.hEC
affinity_cover = create_width_balanced_cover(min(affinity_projection), max(affinity_projection), 4, percent_overlap)

single_edit_affinity_mapptamer = create_mapptamer_graph(edit_dists, affinity_projection, affinity_cover, "single")
single_tree_affinity_mapptamer = create_mapptamer_graph(tree_dists, affinity_projection, affinity_cover, "single")
complete_edit_affinity_mapptamer = create_mapptamer_graph(edit_dists, affinity_projection, affinity_cover, "complete")
complete_tree_affinity_mapptamer = create_mapptamer_graph(tree_dists, affinity_projection, affinity_cover, "complete")

cytoviz1(single_edit_affinity_mapptamer, "1d_affinity_lens_single_linkage_edit_distance")
cytoviz1(single_tree_affinity_mapptamer, "1d_affinity_lens_single_linkage_tree_distance")
cytoviz1(complete_edit_affinity_mapptamer, "1d_affinity_lens_complete_linkage_edit_distance")
cytoviz1(complete_tree_affinity_mapptamer, "1d_affinity_lens_complete_linkage_tree_distance")

selex_projection = node_data$Log2.R3.9
selex_cover = create_width_balanced_cover(min(selex_projection), max(selex_projection), 8, percent_overlap)

single_edit_selex_mapptamer = create_mapptamer_graph(edit_dists, selex_projection, selex_cover, "single")
single_tree_selex_mapptamer = create_mapptamer_graph(tree_dists, selex_projection, selex_cover, "single")
complete_edit_selex_mapptamer = create_mapptamer_graph(edit_dists, selex_projection, selex_cover, "complete")
complete_tree_selex_mapptamer = create_mapptamer_graph(tree_dists, selex_projection, selex_cover, "complete")

cytoviz2(single_edit_selex_mapptamer, "1d_selex_enrichment_lens_single_linkage_edit_distance")
cytoviz2(single_tree_selex_mapptamer, "1d_selex_enrichment_lens_single_linkage_tree_distance")
cytoviz2(complete_edit_selex_mapptamer, "1d_selex_enrichment_lens_complete_linkage_edit_distance")
cytoviz2(complete_tree_selex_mapptamer, "1d_selex_enrichment_lens_complete_linkage_tree_distance")
