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

edit_edges = as.data.frame(read.csv("data/edit_dist_edges.csv"))
tree_edges = as.data.frame(read.csv("data/tree_dist_edges.csv"))
node_data = as.data.frame(read.csv("data/aptamer-data-raw.csv"))

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
