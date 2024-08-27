#' CygnusCluster: Perform Clustering on Dimensionality Reduced Data
#'
#' Clusters data in the `CygnusObject` using specified dimensionality reduction techniques and K-means clustering.
#' The function updates the `CygnusObject` with the resulting cluster assignments.
#'
#' @param CygnusObject A `CygnusObject` containing the data matrices and metadata for analysis.
#' @param use.dims A character string specifying the dimensionality reduction technique to use.
#' Valid options are "PCA" (default), "tSNE", or "UMAP".
#' @param use.markers A logical value indicating whether to use marker genes for clustering (default is TRUE).
#' @param relevant_markers A logical value indicating whether to use only relevant markers for clustering
#' (default is TRUE).
#' @param matrix_name A character string specifying the name of the matrix to use from the `CygnusObject`.
#' The default is "scaled_exp_matrix".
#' @param n_clusters An integer specifying the number of clusters to form (default is 10).
#'
#' @return A `CygnusObject` with the cluster assignments added to the `ev_meta` slot.
#'
#' @details
#' The function performs the following steps:
#' 1. Extracts the specified expression matrix from the `CygnusObject`.
#' 2. Filters the matrix to include only relevant markers if specified.
#' 3. Retrieves the selected dimensionality reduction results (`PCA`, `tSNE`, or `UMAP`).
#' 4. Uses the reduced dimensions for K-means clustering to assign each sample to a cluster.
#' 5. Stores the cluster assignments in the `ev_meta` slot of the `CygnusObject`.
#'
#' @export
ClusterCygnus <- function(CygnusObject,
                          use.dims = "PCA",
                          use.markers = TRUE,
                          relevant_markers = TRUE,
                          matrix_name = "scaled_exp_matrix",
                          n_clusters = 10) {

  expression_matrix <- CygnusObject@matrices[[matrix_name]]

  if (use.markers && relevant_markers) {
    relevant_markers <- CygnusObject@markers_meta$marker[CygnusObject@markers_meta$relevant]
    if (length(relevant_markers) == 0) {
      stop("No relevant markers specified.")
    }
    expression_matrix <- expression_matrix[, relevant_markers, drop = FALSE]
  }

  if (use.dims == "PCA") {
    dim_reduction <- CygnusObject@dim_red$PCA

  } else if (use.dims == "tSNE") {
    dim_reduction <- CygnusObject@dim_red$tSNE

  } else if (use.dims == "UMAP") {
    dim_reduction <- CygnusObject@dim_red$UMAP

  } else {
    stop("Invalid choice for use.dims. Choose 'PCA', 'tSNE', or 'UMAP'.")
  }

  clustering <- kmeans(dim_reduction, centers = n_clusters)$cluster

  CygnusObject@ev_meta$k_means_clusters <- clustering

  return(CygnusObject)
}
