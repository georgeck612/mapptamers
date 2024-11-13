# install mappeR software
require(devtools)
install_github("https://github.com/Uiowa-Applied-Topology/mappeR/tree/cytoscape")
# install.packages("mappeR")
library(mappeR)

# read aptamer data as dataframes
edit_edges = as.data.frame(read.csv("data/editdists2.csv"))
tree_edges = as.data.frame(read.csv("data/treedists2.csv"))
node_data = as.data.frame(read.csv("data/all_aptamers.csv"))

num_datapoints = nrow(node_data)

# create empty edit distance matrix
edit_dists = matrix(0, nrow = num_datapoints, ncol = num_datapoints)
row_indices = match(edit_edges$source, node_data$name)
col_indices = match(edit_edges$target, node_data$name)

# populate edit distance matrix
edit_dists[cbind(row_indices, col_indices)] = edit_edges$editDistance
edit_dists[cbind(col_indices, row_indices)] = edit_edges$editDistance
row.names(edit_dists) = node_data$name
colnames(edit_dists) = node_data$name

# create empty tree distance matrix
tree_dists = matrix(0, nrow = num_datapoints, ncol = num_datapoints)
row_indices = match(tree_edges$source, node_data$name)
col_indices = match(tree_edges$target, node_data$name)

# populate tree distance matrix
tree_dists[cbind(row_indices, col_indices)] = tree_edges$treeDistance
tree_dists[cbind(col_indices, row_indices)] = tree_edges$treeDistance
row.names(tree_dists) = node_data$name
colnames(tree_dists) = node_data$name

row.names(node_data) = node_data$name

ecc = eccentricity_filter(tree_dists)
ecc2 = eccentricity_filter(edit_dists)

# run mapper (clusterball flavor) with varying ball radii
mapptamer1 = create_clusterball_mapper_object(node_data, tree_dists, edit_dists, 40, "single")
mapptamer2 = create_clusterball_mapper_object(node_data, tree_dists, edit_dists, 20, "single")
mapptamer3 = create_clusterball_mapper_object(node_data, tree_dists, edit_dists, 10, "single")
mapptamer4 = create_clusterball_mapper_object(node_data, tree_dists, edit_dists, 5, "single")
mapptamer5 = create_1D_mapper_object(node_data, tree_dists, ecc, create_width_balanced_cover(min(ecc), max(ecc), 5, 10), "single")
mapptamer6 = create_1D_mapper_object(node_data, edit_dists, ecc2, create_width_balanced_cover(min(ecc2), max(ecc2), 5, 10), "single")

# visualize in Cytoscape
cymapper(mapptamer1)
cymapper(mapptamer2)
cymapper(mapptamer3)
cymapper(mapptamer4)
cymapper(mapptamer5)
cymapper(mapptamer6)
