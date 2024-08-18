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
