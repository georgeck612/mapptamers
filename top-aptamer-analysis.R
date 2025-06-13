if(!require("BiocManager")) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install("RCy3", force=TRUE)
}

library(devtools)
install_github("https://github.com/Uiowa-Applied-Topology/mappeR/tree/dev")

library(mappeR)
library(RCy3)

source("draw_trees.R")
source("aptamer_clusterer.R")
source("cytoscape-stuff.R")

# data input and distance matrix creation ---------------------------------

# I don't think there are actually missing values but I can't trust myself
node_data = na.omit(as.data.frame(read.csv("data/all_selex_data.csv")))
edge_data = na.omit(as.data.frame(read.csv("data/aptamer_edges.csv")))

# making sure R doesn't have any cause to yell at me
all_aptamer_names = gsub(">hVSMC-", "Ruiz-RNA_", node_data$Name)
rownames(node_data) = all_aptamer_names
node_data$Name = all_aptamer_names

# flatten interactions column into a list of strings formatted as
# ["source aptamer", "to", "target aptamer"]
squished_edges = unlist(strsplit(edge_data$name, " "))

# recover sources/targets via above formatting
sources = squished_edges[seq(1, length(squished_edges), 3)]
targets = squished_edges[seq(3, length(squished_edges), 3)]

edge_data$sources = sources
edge_data$targets = targets

missing_aptamers = setdiff(all_aptamer_names, unique(union(sources, targets)))

our_aptamers = setdiff(all_aptamer_names, missing_aptamers)

print(nrow(edge_data))
edge_data = edge_data[edge_data$source %in% our_aptamers,]
edge_data = edge_data[edge_data$target %in% our_aptamers,]
print(nrow(edge_data))

print(nrow(node_data))
node_data = node_data[node_data$Name %in% our_aptamers,]
print(nrow(node_data))

sources = edge_data$sources
targets = edge_data$targets

num_aptamers = length(unique(union(sources, targets)))

# grab entire source/target aptamer locations in mother dataset
row_indices = match(sources, node_data$Name)
col_indices = match(targets, node_data$Name)

# create empty edit distance matrix
edit_dists = matrix(0, nrow = num_aptamers, ncol = num_aptamers)

# populate edit distance matrix
edit_dists[cbind(row_indices, col_indices)] = edge_data$editDistance
edit_dists[cbind(col_indices, row_indices)] = edge_data$editDistance

# again I'm being paranoid
rownames(edit_dists) = node_data$Name
colnames(edit_dists) = node_data$Name

# create empty tree distance matrix
tree_dists = matrix(0, nrow = num_aptamers, ncol = num_aptamers)

# populate tree distance matrix
tree_dists[cbind(row_indices, col_indices)] = edge_data$treeDistance
tree_dists[cbind(col_indices, row_indices)] = edge_data$treeDistance

# again I'm being paranoid
rownames(tree_dists) = node_data$Name
colnames(tree_dists) = node_data$Name


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
create_mapptamer_graph <- function(dists, filtered, cover, linkage_method, metric) {
  # mapper time
  mapptamer = create_1D_mapper_object(node_data, dists, filtered, cover, clusterer = global_hierarchical_clusterer(linkage_method, dists, metric))

  # get aptamers in vertices of mapper graph
  aptamer_balls = lapply(mapptamer[[1]]$data, function(x) node_data[unlist(strsplit(x, ", ")),])

  # calculate median selex data across mapper vertices
  median_log10_selex_enrichment = sapply(aptamer_balls, function(x) median(log10(node_data[x$Name, "RP10M.9"]/node_data[x$Name, "RP10M.3A"])))

  # calculate median ASSET data across mapper vertices
  median_log2_human_affinity = sapply(aptamer_balls, function(x) median(log2(node_data[x$Name, "ASSET.hVSMC.hEC"])))
  median_log2_mouse_affinity = sapply(aptamer_balls, function(x) median(log2(node_data[x$Name, "ASSET.mVSMC.mEC"])))

  # attach calculated info to mapper dataframe
  mapptamer[[1]]$median_log10_selex_enrichment = median_log10_selex_enrichment
  mapptamer[[1]]$median_log2_human_affinity = median_log2_human_affinity
  mapptamer[[1]]$median_log2_mouse_affinity = median_log2_mouse_affinity

  return(mapptamer)
}


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

# NODE SIZE: selex score
# NODE BORDER COLOR: human affinity
# NODE FILL COLOR: average distance to medoid
cytoviz1 = function(mapptamer, name) {
  size_col = "median_log10_selex_enrichment"
  size_data = sort(mapptamer[[1]][, size_col])
  border_color_col = "median_log2_human_affinity"
  border_color_data = sort(mapptamer[[1]][, border_color_col])
  fill_color_col = "mean_dist_to_medoid"
  fill_color_data = sort(mapptamer[[1]][, fill_color_col])
  visualize_mapper_data(mapptamer, size_col, size_data, border_color_col, border_color_data, fill_color_col, fill_color_data, name)
}

# NODE SIZE: human affinity
# NODE BORDER COLOR: selex score
# NODE FILL COLOR: average distance to medoid
cytoviz2 = function(mapptamer, name) {
  size_col = "median_log2_human_affinity"
  size_data = sort(mapptamer[[1]][, size_col])
  border_color_col = "median_log10_selex_enrichment"
  border_color_data = sort(mapptamer[[1]][, border_color_col])
  fill_color_col = "mean_dist_to_medoid"
  fill_color_data = sort(mapptamer[[1]][, fill_color_col])
  visualize_mapper_data(mapptamer, size_col, size_data, border_color_col, border_color_data, fill_color_col, fill_color_data, name)
}

percent_overlap = 25
affinity_projection = log2(node_data$ASSET.hVSMC.hEC)
affinity_cover = create_width_balanced_cover(min(affinity_projection), max(affinity_projection), 8, percent_overlap)

single_edit_affinity_mapptamer = create_mapptamer_graph(edit_dists, affinity_projection, affinity_cover, "single", "Edit")
single_tree_affinity_mapptamer = create_mapptamer_graph(tree_dists, affinity_projection, affinity_cover, "single", "Tree")

cytoviz1(single_edit_affinity_mapptamer, "1d_affinity_lens_single_linkage_edit_distance")
cytoviz1(single_tree_affinity_mapptamer, "1d_affinity_lens_single_linkage_tree_distance")

selex_projection = log10(node_data$RP10M.9/node_data$RP10M.3A)
selex_cover = create_width_balanced_cover(min(selex_projection), max(selex_projection), 10, percent_overlap)

single_edit_selex_mapptamer = create_mapptamer_graph(edit_dists, selex_projection, selex_cover, "single", "Edit")
single_tree_selex_mapptamer = create_mapptamer_graph(tree_dists, selex_projection, selex_cover, "single", "Tree")

cytoviz2(single_edit_selex_mapptamer, "1d_selex_enrichment_lens_single_linkage_edit_distance")
cytoviz2(single_tree_selex_mapptamer, "1d_selex_enrichment_lens_single_linkage_tree_distance")
