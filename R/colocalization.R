#' Generate an UpSet Plot for Relevant Markers
#'
#' @param data An object of class \code{CygnusObject} containing the binary expression matrix.
#' @param markers Character vector specifying which markers to include in the UpSet plot. If NULL, all relevant markers are used.
#' @param ... Additional arguments to pass to the upset function.
#' @return An UpSet plot visualizing the intersections of relevant markers.
#' @export
plotUpSet <- function(data, markers = NULL, ...) {
  # Check if binary_exp_matrix exists
  if (!"binary_exp_matrix" %in% names(data@matrices)) {
    stop("Binary expression matrix not found. Please run createBinaryMatrix first.")
  }

  binary_matrix <- data@matrices$binary_exp_matrix

  relevant_markers <- data@markers_meta$marker[data@markers_meta$relevant]
  if (is.null(relevant_markers) || length(relevant_markers) == 0) {
    stop("No relevant markers found.")
  }

  if (!is.null(markers)) {
    relevant_markers <- relevant_markers[relevant_markers %in% markers]
    if (length(relevant_markers) == 0) {
      stop("None of the specified markers were found in relevant markers.")
    }
  }

  binary_matrix <- binary_matrix[, relevant_markers, drop = FALSE]

  binary_df <- as.data.frame(binary_matrix)
  colnames(binary_df) <- relevant_markers

  UpSetR::upset(binary_df, ...)
}
