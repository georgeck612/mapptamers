# visualize_mapper_data <- function(mapper_data, is_ballmapper = TRUE) {
#   nodes = mapper_data[[1]]
#   edges = mapper_data[[2]]
#
#   createNetworkFromDataFrames(nodes, edges)
#
#   style.name = paste("mapperstyle", runif(1))
#   defaults <- list(
#     NODE_SHAPE = "ellipse",
#     NODE_BORDER_WIDTH = 10,
#     NODE_BORDER_PAINT = "#000"
#   )
#
#   nodeSizes <- mapVisualProperty('node size',
#                                  'id',
#                                  'd',
#                                  1:nrow(nodes),
#                                  100 * sqrt(nodes$cluster_size / max(nodes$cluster_size)))
#   edgeWidth <- mapVisualProperty('edge width', 'weight', 'c', c(0, .5, 1), c(0, 5, 10))
#
#   fill_colors = lapply(nodes$tightness / max(nodes$tightness), function(x)
#     rgb(x, x, x))
#
#   nodeFillColors <- mapVisualProperty('node fill color', 'id', 'd', 1:nrow(nodes), fill_colors)
#
#   if (is_ballmapper) {
#     # ballmapper needs no more styling
#     createVisualStyle(style.name,
#                       defaults,
#                       list(nodeSizes, edgeWidth, nodeFillColors))
#   } else {
#     # conventional mapper needs bin coloring
#     num_bins = length(unique(nodes$bin))
#     bin_colors = viridisLite::plasma(num_bins)
#     nodeBorderColors <- mapVisualProperty('node border color', 'bin', 'd', 1:num_bins, bin_colors)
#     createVisualStyle(
#       style.name,
#       defaults,
#       list(nodeSizes, edgeWidth, nodeBorderColors, nodeFillColors)
#     )
#   }
#
#   setVisualStyle(style.name)
# }

# Cytoscape ---------------------------------------------------------------

library(viridisLite)

visualize_mapper_data <- function(mapper_data, size_col, size_data, border_color_col, border_color_data, fill_color_col, fill_color_data, name) {
  nodes = mapper_data[[1]]
  edges = mapper_data[[2]]

  createNetworkFromDataFrames(nodes, edges, title = name)

  style.name = paste("selexstyle", runif(1))

  defaults <- list(
    NODE_SHAPE = "ellipse",
    NODE_BORDER_WIDTH = 25
  )

  nodeSizes <- mapVisualProperty('node size', size_col, 'd', size_data, ((7*size_data / max(size_data)) + 5)^2)
  edgeWidth <- mapVisualProperty('edge width', 'weight', 'c', c(0, .5, 1), c(0, 5, 10))
  nodeFillColors <- mapVisualProperty('node fill color', fill_color_col, 'c', c(min(fill_color_data), max(fill_color_data)), c('black', 'white'))
  nodeBorderColors <- mapVisualProperty('node border paint', border_color_col, 'd', border_color_data, lapply(plasma(length(border_color_data)), function(x) substr(x, 1, nchar(x)-2)))
  createVisualStyle(
    style.name,
    defaults,
    list(nodeSizes, edgeWidth, nodeBorderColors, nodeFillColors)
  )
  setVisualStyle(style.name)
  layoutNetwork('attributes-layout nodeAttribute=patch')
}



#' Open mapper graph in Cytoscape
#'
#' @param mapperobject A set of data frames representing a mapper object, returned by, say, [create_mapper_object()].
#'
#' @return Nothing; opens Cytoscape with information from the mapper object ported there. Cytoscape must be actively running in the background for this method to work.
#' @export
#'
#' @examples
#' \dontrun{
#' # this example requires Cytoscape to be open and running in the background to work properly
#'
#' data = data.frame(x = sapply(1:100, function(x) cos(x)), y = sapply(1:100, function(x) sin(x)))
#' projx = data$x
#'
#' num_bins = 10
#' percent_overlap = 25
#'
#' cover = get_width_balanced_cover(min(projx), max(projx), num_bins, percent_overlap)
#'
#' mapperobj = create_1D_mapper_object(data, dist(data), projx, cover, "single")
#' cymapper(mapperobj)
#' }
cymapper <- function(mapper_data, size_attr, border_color_attr, fill_color_attr) {
  # pass to visualizer for........visualizing...
  visualize_mapper_data(mapper_data, size_col, size_data, border_color_col, border_color_data, fill_color_col, fill_color_data)

  # if this isn't here R will print something useless
  return(invisible(NULL))
}
