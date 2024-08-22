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
    data@markers_meta$marker <- colnames(test@matrices$Raw_Score)
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
#' @param seed Numeric value for setting a random seed. Default is 42.
#' @param ... Additional arguments to pass to the Rtsne function.
#' @return The updated CygnusObject with t-SNE coordinates stored in the 'dim_red' slot.
#' @export
runTSNE <- function(data, matrix_name = "Raw_Score", use_relevant = TRUE, dims = 2, use_pcs = FALSE, seed = 42, ...) {
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

  # Set the random seed for reproducibility
  set.seed(seed)

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
#' @param n_neighbors Numeric value specifying the number of neighbors for UMAP. Default is 30.
#' @param min_dist Numeric value specifying the minimum distance for UMAP. Default is 0.3.
#' @param verbose Logical value indicating whether to display a progress bar. Default is TRUE.
#' @param ... Additional arguments to pass to the \code{uwot::umap} function.
#' @return The updated CygnusObject with UMAP coordinates stored in the 'dim_red' slot.
#' @export
runUMAP <- function(data, matrix_name = "Raw_Score", use_relevant = TRUE, dims = 3, use_pcs = FALSE,
                    n_neighbors = 30, min_dist = 0.3, verbose = TRUE, ...) {
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
      uwot::umap(expression_matrix, n_neighbors = n_neighbors, min_dist = min_dist,
                 n_components = dims, verbose = verbose, ...)
    },
    error = function(e) {
      stop("Error in UMAP computation: ", e$message)
    }
  )

  data@dim_red$UMAP <- umap_result

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
#' @param ... Additional arguments for the plot function
#' @return A plot of the PCA results.
#' @export
plotPCA <- function(data, plot_3d = FALSE, color_by = NULL, ...) {
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

    color_palette <- RColorBrewer::brewer.pal(min(length(levels(plot_data$meta)), 9), "Set1")
    plot_data$color <- color_palette[plot_data$meta]

    legend_colors <- setNames(color_palette[1:length(levels(plot_data$meta))], levels(plot_data$meta))
  } else if (!is.null(color_by)) {
    warning(paste("Metadata column", color_by, "not found. Plotting without color."))
  }

  if (plot_3d && ncol(pca_coords) >= 3) {
    if (!is.null(plot_data$color)) {
      plotly::plot_ly(plot_data, x = ~V1, y = ~V2, z = ~V3, color = ~meta, colors = color_palette, type = 'scatter3d', mode = 'markers') %>%
        plotly::layout(legend = list(title = list(text = color_by)))
    } else {
      plotly::plot_ly(plot_data, x = ~V1, y = ~V2, z = ~V3, type = 'scatter3d', mode = 'markers')
    }
  } else {
    # Set up the layout to leave space on the right for the legend
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))

    par(mar = c(5, 4, 4, 10) + 0.1)  # Increase right margin for legend space

    if (!is.null(plot_data$color)) {
      plot(plot_data$V1, plot_data$V2, col = plot_data$color, xlab = "PC1", ylab = "PC2", main = "PCA", ...)
      legend("topright", inset = c(-0.2, 0), legend = names(legend_colors), fill = legend_colors, title = color_by, xpd = TRUE)
    } else {
      plot(plot_data$V1, plot_data$V2, xlab = "PC1", ylab = "PC2", main = "PCA", ...)
    }
  }
}


#' Visualize UMAP Results
#'
#' @param data An object of class \code{CygnusObject} containing UMAP results.
#' @param plot_3d Logical value indicating whether to plot in 3D. Default is FALSE.
#' @param color_by Character string specifying metadata column to color by. Default is NULL (no coloring).
#' @param marker_size Numeric value specifying the size of the markers in the plot. Default is 1.
#' @return A plot of the UMAP results.
#' @export
plotUMAP <- function(data, plot_3d = FALSE, color_by = NULL, marker_size = 1) {
  if (!"UMAP" %in% names(data@dim_red)) {
    stop("UMAP results not found. Run UMAP first.")
  }

  umap_coords <- data@dim_red$UMAP
  if (is.null(umap_coords) || nrow(umap_coords) == 0 || ncol(umap_coords) < 2) {
    stop("UMAP coordinates are not available or insufficient for plotting.")
  }

  plot_data <- as.data.frame(umap_coords)
  colnames(plot_data) <- paste0("V", 1:ncol(plot_data))

  if (!is.null(color_by) && color_by %in% names(data@ev_meta)) {
    plot_data$meta <- data@ev_meta[[color_by]]
    plot_data$meta <- as.factor(plot_data$meta)

    if (is.factor(plot_data$meta)) {
      color_palette <- RColorBrewer::brewer.pal(length(levels(plot_data$meta)), "Set1")
      plot_data$color <- color_palette[plot_data$meta]
    }
  } else if (!is.null(color_by)) {
    warning(paste("Metadata column", color_by, "not found. Plotting without color."))
  }

  if (plot_3d && ncol(umap_coords) >= 3) {
    if (!is.null(plot_data$color)) {
      plotly::plot_ly(plot_data, x = ~V1, y = ~V2, z = ~V3, color = ~meta, type = 'scatter3d', mode = 'markers', marker = list(size = marker_size))
    } else {
      plotly::plot_ly(plot_data, x = ~V1, y = ~V2, z = ~V3, type = 'scatter3d', mode = 'markers', marker = list(size = marker_size))
    }
  } else {
    if (!is.null(plot_data$color)) {
      plot(plot_data$V1, plot_data$V2, col = plot_data$color, xlab = "UMAP1", ylab = "UMAP2", main = "UMAP", pch = 16, cex = marker_size)
    } else {
      plot(plot_data$V1, plot_data$V2, xlab = "UMAP1", ylab = "UMAP2", main = "UMAP", pch = 16, cex = marker_size)
    }
  }
}

#' Visualize t-SNE Results
#'
#' @param data An object of class \code{CygnusObject} containing t-SNE results.
#' @param plot_3d Logical value indicating whether to plot in 3D. Default is FALSE.
#' @param color_by Character string specifying metadata column to color by. Default is NULL (no coloring).
#' @param marker_size Numeric value specifying the size of the markers in the plot. Default is 1.
#' @return A plot of the t-SNE results.
#' @export
plotTSNE <- function(data, plot_3d = FALSE, color_by = NULL, marker_size = 1, ...) {
  if (!"tSNE" %in% names(data@dim_red)) {
    stop("t-SNE results not found. Run t-SNE first.")
  }

  tsne_coords <- data@dim_red$tSNE
  if (is.null(tsne_coords) || nrow(tsne_coords) == 0 || ncol(tsne_coords) < 2) {
    stop("t-SNE coordinates are not available or insufficient for plotting.")
  }

  plot_data <- as.data.frame(tsne_coords)
  colnames(plot_data) <- paste0("V", 1:ncol(plot_data))

  if (!is.null(color_by) && color_by %in% names(data@ev_meta)) {
    plot_data$meta <- data@ev_meta[[color_by]]
    plot_data$meta <- as.factor(plot_data$meta)

    if (is.factor(plot_data$meta)) {
      color_palette <- RColorBrewer::brewer.pal(length(levels(plot_data$meta)), "Set1")
      plot_data$color <- color_palette[plot_data$meta]
    }
  } else if (!is.null(color_by)) {
    warning(paste("Metadata column", color_by, "not found. Plotting without color."))
  }

  if (plot_3d && ncol(tsne_coords) >= 3) {
    if (!is.null(plot_data$color)) {
      plotly::plot_ly(plot_data, x = ~V1, y = ~V2, z = ~V3, color = ~meta, type = 'scatter3d', mode = 'markers', marker = list(size = marker_size)) %>%
        plotly::layout(legend = list(title = list(text = color_by)))
    } else {
      plotly::plot_ly(plot_data, x = ~V1, y = ~V2, z = ~V3, type = 'scatter3d', mode = 'markers', marker = list(size = marker_size))
    }
  } else {
    if (!is.null(plot_data$color)) {
      plot(plot_data$V1, plot_data$V2, col = plot_data$color, xlab = "t-SNE1", ylab = "t-SNE2", main = "t-SNE", pch = 16, cex = marker_size, ...)
      legend("topright", legend = levels(plot_data$meta), fill = unique(plot_data$color), title = color_by, cex = 0.8)
    } else {
      plot(plot_data$V1, plot_data$V2, xlab = "t-SNE1", ylab = "t-SNE2", main = "t-SNE", pch = 16, cex = marker_size, ...)
    }
  }
}


