#' Plot Distributions of Selected Markers
#'
#' This function plots the distribution of expression levels for selected markers in an expression matrix.
#' The user can specify a subset of markers to plot, or plot all markers by default.
#'
#' @param obj An object containing the expression matrix, typically a Seurat or similar object.
#' @param plot_markers A vector of marker names to plot. If set to "ALL", distributions for all markers will be plotted. Default is "ALL".
#' @param matrix Character string specifying the name of the slot containing the expression matrix in \code{obj}. Default is "exp_matrix".
#' @return A series of histogram plots, each representing the distribution of expression levels for the specified markers.
#' @export
plotDistribution <- function(
    obj,
    plot_markers = "ALL",
    matrix = "exp_matrix"
) {
  if (plot_markers == "ALL") {
    plot_markers <- colnames(obj@exp_matrix)
  } else {
    plot_markers <- intersect(plot_markers, colnames(obj@exp_matrix))
  }

  num_cols <- length(plot_markers)
  num_rows <- ceiling(num_cols / 3)

  par(mfrow = c(num_rows, 3))
  par(mar = c(2, 2, 2, 1))

  for (marker in plot_markers) {
    hist(obj@exp_matrix[, marker], breaks = 1000,
         main = paste(marker),
         xlab = "", ylab = "")
  }

  par(mfrow = c(1, 1))
}

#' Plot Average Expression Heatmap by Group
#'
#' This function generates a heatmap of average marker expressions, grouped by a specified metadata column.
#' The heatmap allows customization of clustering distance, color palette, font size, and scaling method.
#'
#' @param data An object containing an expression matrix and metadata, typically a Seurat or similar object.
#' @param group_column Character string specifying the metadata column to group by.
#' @param clustering_distance Character string specifying the distance metric for clustering columns. Default is "euclidean".
#' @param colors A color palette to use for the heatmap. Default is a reversed "RdYlBu" palette from RColorBrewer.
#' @param fontsize Numeric value specifying the font size for heatmap text. Default is 8.
#' @param scale Character string specifying whether the data should be scaled by 'row', 'column', or 'none'. Default is 'column'.
#' @param cluster_rows Logical value indicating whether to cluster rows in the heatmap. Default is FALSE.
#' @param na.rm Logical value indicating whether to remove NA values before calculating averages. Default is TRUE.
#' @return A heatmap plot showing average marker expressions grouped by the specified metadata column.
#' @export
plotAvgHeatmap <- function(data, group_column,
                           clustering_distance = "euclidean",
                           colors = rev(colorRampPalette(RColorBrewer::brewer.pal(n = 11, "RdYlBu"))(100)),
                           fontsize = 8,
                           scale = 'column',
                           cluster_rows = FALSE,
                           na.rm = TRUE) {

  expression_matrix <- data@exp_matrix
  group_metadata <- data@metadata[[group_column]]

  combined_data <- cbind(expression_matrix, group_metadata)

  avg_marker_expressions <- combined_data %>%
    group_by(group_metadata) %>%
    summarise(across(everything(), mean, na.rm = na.rm)) %>%
    select(-group_metadata) %>%
    as.matrix()

  pheatmap(avg_marker_expressions,
           clustering_distance_cols = clustering_distance,
           cluster_rows = cluster_rows,
           color = colors,
           main = paste("Average Marker Expressions by", group_column),
           fontsize = fontsize,
           scale = scale,
           labels_row = unique(group_metadata))
}

#' Scale Expression Matrix by Maximum Value
#'
#' This function scales the expression matrix by dividing each marker's values by the maximum value in that marker.
#' The scaled matrix is added as an additional layer in the data object.
#'
#' @param data An object containing an expression matrix and metadata, typically a Seurat or similar object.
#' @param matrix_name Character string specifying the name for the new scaled matrix layer. Default is "scaled_exp_matrix".
#' @return The data object with an additional layer containing the scaled expression matrix.
#' @export
scaleExpressionMatrix <- function(data, matrix_name = "scaled_exp_matrix") {

  expression_matrix <- data@exp_matrix

  max_values <- apply(expression_matrix, 2, max, na.rm = TRUE)
  scaled_matrix <- sweep(expression_matrix, 2, max_values, `/`)

  data[[matrix_name]] <- scaled_matrix

  return(data)
}


