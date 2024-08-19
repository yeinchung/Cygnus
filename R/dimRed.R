#' Perform PCA on a Specified Matrix and Store in Dimensionality Reduction Slot
#'
#' @param data An object of class \code{CygnusObject} containing an expression matrix.
#' @param matrix_name Character string specifying the name of the matrix to use. Default is "Raw_Score".
#' @param num_components Numeric value specifying the number of principal components to compute. Default is 50.
#' @param use_relevant Logical indicating whether to use only relevant markers. Default is TRUE.
#' @return The updated CygnusObject with PCA coordinates stored in the 'dim_red' slot.
#' @export
runPCA <- function(data, matrix_name = "Raw_Score", num_components = 50, use_relevant = TRUE) {
  if (!(matrix_name %in% names(data@matrices))) {
    stop(paste("Matrix", matrix_name, "not found in CygnusObject"))
  }

  expression_matrix <- data@matrices[[matrix_name]]

  if (use_relevant) {
    relevant_markers <- data@markers_meta$marker[data@markers_meta$relevant]
    if (length(relevant_markers) == 0) {
      stop("No relevant markers specified.")
    }
    expression_matrix <- expression_matrix[, relevant_markers, drop = FALSE]
  }

  if (nrow(expression_matrix) == 0 || ncol(expression_matrix) == 0) {
    stop("The expression matrix has zero rows or columns.")
  }

  pca_result <- tryCatch(
    {
      prcomp(expression_matrix, center = TRUE, scale. = TRUE)
    },
    error = function(e) {
      stop("Error in PCA computation: ", e$message)
    }
  )

  available_components <- min(num_components, ncol(pca_result$x))
  data@dim_red$PCA <- pca_result$x[, 1:available_components, drop = FALSE]

  return(data)
}

#' Perform t-SNE on a Specified Matrix and Store in Dimensionality Reduction Slot
#'
#' @param data An object of class \code{CygnusObject} containing an expression matrix.
#' @param matrix_name Character string specifying the name of the matrix to use. Default is "Raw_Score".
#' @param use_relevant Logical indicating whether to use only relevant markers. Default is TRUE.
#' @param dims Numeric value specifying the number of dimensions for t-SNE. Default is 2.
#' @param use_pcs Logical indicating whether to use PCA results as input. Default is FALSE.
#' @param ... Additional arguments to pass to the Rtsne function.
#' @return The updated CygnusObject with t-SNE coordinates stored in the 'dim_red' slot.
#' @export
runTSNE <- function(data, matrix_name = "Raw_Score", use_relevant = TRUE, dims = 2, use_pcs = FALSE, ...) {
  if (!(matrix_name %in% names(data@matrices))) {
    stop(paste("Matrix", matrix_name, "not found in CygnusObject"))
  }

  if (use_pcs && !"PCA" %in% names(data@dim_red)) {
    stop("PCA results not found. Run PCA first.")
  }

  expression_matrix <- data@matrices[[matrix_name]]

  if (use_relevant) {
    relevant_markers <- data@markers_meta$marker[data@markers_meta$relevant]
    if (length(relevant_markers) == 0) {
      stop("No relevant markers specified.")
    }
    expression_matrix <- expression_matrix[, relevant_markers, drop = FALSE]
  }

  if (use_pcs) {
    expression_matrix <- data@dim_red$PCA
  }

  if (nrow(expression_matrix) == 0 || ncol(expression_matrix) == 0) {
    stop("The expression matrix has zero rows or columns.")
  }

  # need to add random seed
  tsne_result <- tryCatch(
    {
      Rtsne::Rtsne(expression_matrix, dims = dims, ...)
    },
    error = function(e) {
      stop("Error in t-SNE computation: ", e$message)
    }
  )

  data@dim_red$tSNE <- tsne_result$Y

  return(data)
}

#' Perform UMAP on a Specified Matrix and Store in Dimensionality Reduction Slot
#'
#' @param data An object of class \code{CygnusObject} containing an expression matrix.
#' @param matrix_name Character string specifying the name of the matrix to use. Default is "Raw_Score".
#' @param use_relevant Logical indicating whether to use only relevant markers. Default is TRUE.
#' @param dims Numeric value specifying the number of dimensions for UMAP. Default is 2.
#' @param use_pcs Logical indicating whether to use PCA results as input. Default is FALSE.
#' @param ... Additional arguments to pass to the umap function.
#' @return The updated CygnusObject with UMAP coordinates stored in the 'dim_red' slot.
#' @export
runUMAP <- function(data, matrix_name = "Raw_Score", use_relevant = TRUE, dims = 2, use_pcs = FALSE, ...) {
  if (!(matrix_name %in% names(data@matrices))) {
    stop(paste("Matrix", matrix_name, "not found in CygnusObject"))
  }

  if (use_pcs && !"PCA" %in% names(data@dim_red)) {
    stop("PCA results not found. Run PCA first.")
  }

  expression_matrix <- data@matrices[[matrix_name]]

  if (use_relevant) {
    relevant_markers <- data@markers_meta$marker[data@markers_meta$relevant]
    if (length(relevant_markers) == 0) {
      stop("No relevant markers specified.")
    }
    expression_matrix <- expression_matrix[, relevant_markers, drop = FALSE]
  }

  if (use_pcs) {
    expression_matrix <- data@dim_red$PCA
  }


  if (nrow(expression_matrix) == 0 || ncol(expression_matrix) == 0) {
    stop("The expression matrix has zero rows or columns.")
  }


  umap_result <- tryCatch(
    {
      umap::umap(expression_matrix, n_neighbors = 15, min_dist = 0.3, n_components = dims, ...)
    },
    error = function(e) {
      stop("Error in UMAP computation: ", e$message)
    }
  )

  data@dim_red$UMAP <- umap_result$layout

  return(data)
}

#' Plot Elbow Plot for PCA
#'
#' @param data An object of class \code{CygnusObject} containing PCA results.
#' @param ... Additional arguments to pass to the plot function.
#' @return An elbow plot showing the variance explained by each principal component.
#' @export
plotElbowPlot <- function(data, ...) {
  if (!"PCA" %in% names(data@dim_red)) {
    stop("PCA results not found. Run PCA first.")
  }

  pca_coordinates <- data@dim_red$PCA
  pca_variance <- apply(pca_coordinates, 2, var)

  plot(pca_variance, type = "b", pch = 19, xlab = "Principal Component", ylab = "Variance Explained",
       main = "Elbow Plot of Principal Components", ...)
}

#' Visualize PCA Results
#'
#' @param data An object of class \code{CygnusObject} containing PCA results.
#' @param plot_3d Logical value indicating whether to plot in 3D. Default is FALSE.
#' @param color_by Character string specifying metadata column to color by. Default is NULL (no coloring).
#' @return A plot of the PCA results.
#' @export
plotPCA <- function(data, plot_3d = FALSE, color_by = NULL) {
  if (!"PCA" %in% names(data@dim_red)) {
    stop("PCA results not found. Run PCA first.")
  }

  pca_coords <- data@dim_red$PCA
  if (is.null(pca_coords) || nrow(pca_coords) == 0 || ncol(pca_coords) < 2) {
    stop("PCA coordinates are not available or insufficient for plotting.")
  }

  plot_data <- as.data.frame(pca_coords)
  colnames(plot_data) <- paste0("V", 1:ncol(plot_data))

  if (!is.null(color_by) && color_by %in% names(data@ev_meta)) {
    plot_data$meta <- data@ev_meta[[color_by]]
    plot_data$meta <- as.factor(plot_data$meta)

    if (is.factor(plot_data$meta)) { # need to fix
      color_palette <- RColorBrewer::brewer.pal(length(levels(plot_data$meta)), "Set1")
      plot_data$color <- color_palette[plot_data$meta]
    }
  } else if (!is.null(color_by)) {
    warning(paste("Metadata column", color_by, "not found. Plotting without color."))
  }

  if (plot_3d && ncol(pca_coords) >= 3) {
    if (!is.null(plot_data$color)) {
      plotly::plot_ly(plot_data, x = ~V1, y = ~V2, z = ~V3, color = ~meta, type = 'scatter3d', mode = 'markers')
    } else {
      plotly::plot_ly(plot_data, x = ~V1, y = ~V2, z = ~V3, type = 'scatter3d', mode = 'markers')
    }
  } else {
    if (!is.null(plot_data$color)) {
      plot(plot_data$V1, plot_data$V2, col = plot_data$color, xlab = "PC1", ylab = "PC2", main = "PCA")
    } else {
      plot(plot_data$V1, plot_data$V2, xlab = "PC1", ylab = "PC2", main = "PCA")
    }
  }
}


