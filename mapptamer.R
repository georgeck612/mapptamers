require(devtools)
install_github("https://github.com/Uiowa-Applied-Topology/mappeR/tree/dev")
library(mappeR)

edit_edges = as.data.frame(read.csv("data/edit_dist_edges.csv"))
tree_edges = as.data.frame(read.csv("data/tree_dist_edges.csv"))
node_data = as.data.frame(read.csv("data/aptamer-data-raw.csv"))

num_datapoints = nrow(node_data)

edit_dists = matrix(0, nrow = num_datapoints, ncol = num_datapoints)
row_indices = match(edit_edges$source, node_data$name)
col_indices = match(edit_edges$target, node_data$name)

edit_dists[cbind(row_indices, col_indices)] = edit_edges$editDistance
edit_dists[cbind(col_indices, row_indices)] = edit_edges$editDistance
row.names(edit_dists) = node_data$name
colnames(edit_dists) = node_data$name

tree_dists = matrix(0, nrow = num_datapoints, ncol = num_datapoints)
row_indices = match(tree_edges$source, node_data$name)
col_indices = match(tree_edges$target, node_data$name)

tree_dists[cbind(row_indices, col_indices)] = tree_edges$treeDistance
tree_dists[cbind(col_indices, row_indices)] = tree_edges$treeDistance
row.names(tree_dists) = node_data$name
colnames(tree_dists) = node_data$name

row.names(node_data) = node_data$name

mapptamer1 = create_clusterball_mapper_object(node_data, tree_dists, edit_dists, 40, "single")
mapptamer2 = create_clusterball_mapper_object(node_data, tree_dists, edit_dists, 20, "single")
mapptamer3 = create_clusterball_mapper_object(node_data, tree_dists, edit_dists, 10, "single")
mapptamer4 = create_clusterball_mapper_object(node_data, tree_dists, edit_dists, 5, "single")

cymapper(mapptamer1)
cymapper(mapptamer2)
cymapper(mapptamer3)
cymapper(mapptamer4)
