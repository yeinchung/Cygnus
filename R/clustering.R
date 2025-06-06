#' Perform Clustering on Dimensionality Reduced Data
#'
#' Clusters data in the `CygnusObject` using specified dimensionality reduction techniques and K-means clustering.
#' The function updates the `CygnusObject` with the resulting cluster assignments.
#'
#' @param CygnusObject A `CygnusObject` containing the data matrices and metadata for analysis.
#' @param use.dims A character string specifying the dimensionality reduction technique to use.
#' Valid options are "PCA" (default), "tSNE", "UMAP", or None.
#' @param relevant_markers A logical value indicating whether to use only relevant markers for clustering
#' (default is TRUE).
#' @param matrix_name A character string specifying the name of the matrix to use from the `CygnusObject`.
#' The default is "scaled_exp_matrix".
#' @param clustering_method A character string specifying the method for clustering.
#' "kmeans", "hdbscan", and "leiden" are available.
#' @param n_clusters An integer specifying the number of clusters to form (default is 3).
#' Only used when clustering_method is set to "kmeans".
#' @param hdb_min_cluster_size An integer specifying the smallest size grouping to be considered as a cluster (default is 10).
#' Only used when clustering_method is set to "hdmscan".
#' @param graph_distance The maximum number of nearest neighbours to compute. Used for leiden.
#' @param leiden_resolution Resolution parameter for leiden.
#'
#' @return A `CygnusObject` with the cluster assignments added to the `ev_meta` slot.
#'
#' @details
#' The function performs the following steps:
#' 1. Extracts the specified expression matrix from the `CygnusObject`.
#' 2. Filters the matrix to include only relevant markers if specified.
#' 3. Retrieves the selected dimensionality reduction results (`PCA`, `tSNE`, or `UMAP` or `None` ).
#' 4. Uses the reduced dimensions for K-means clustering to assign each sample to a cluster.
#' 5. Stores the cluster assignments in the `ev_meta` slot of the `CygnusObject`.
#'
#' @export
# K-means, Leiden , HDBscan


ClusterCygnus <- function(CygnusObject,
                          use.dims = "PCA",
                          relevant_markers = TRUE,
                          matrix_name = "scaled_exp_matrix",
                          clustering_method = "kmeans",
                          n_clusters = 3,
                          hdb_min_cluster_size = 10,
                          graph_distance = 15,
                          leiden_resolution =1.0) {


  # --- Dim reduction
  if (use.dims == "PCA") {
    # clustering_input <- runPCA(CygnusObject, matrix_name = matrix_name)
    clustering_input <- CygnusObject@dim_red$PCA
  }

  else if (use.dims == "tSNE") {
    # clustering_input <- runTSNE(CygnusObject, matrix_name = matrix_name)
    clustering_input <- CygnusObject@dim_red$tSNE
  }

  else if (use.dims == "UMAP") {
    # clustering_input <- runUMAP(CygnusObject, matrix_name = matrix_name)
    clustering_input <- CygnusObject@dim_red$UMAP
  }

  else if (use.dims == "None") {
    clustering_input <- CygnusObject@matrices[[matrix_name]]

    if (relevant_markers) {
      relevant_markers <- CygnusObject@markers_meta$marker[CygnusObject@markers_meta$relevant]
      if (length(relevant_markers) == 0) {
        stop("No relevant markers specified.")
      }
      clustering_input <- clustering_input[, relevant_markers, drop = FALSE]
    }

  }

  else {
    stop("Invalid choice for use.dims. Choose 'PCA', 'tSNE', 'UMAP' or 'None'.")
  }

  if (clustering_method == "kmeans") {
    result_clusters <- kmeans(clustering_input, centers = n_clusters)$cluster
    CygnusObject@ev_meta[['k_means_clusters']] <- result_clusters

  } else if (clustering_method == "hdbscan") {
    hdb_result <- dbscan::hdbscan(clustering_input, minPts = hdb_min_cluster_size)
    CygnusObject@ev_meta[['hdbscan_clusters']] <- hdb_result$cluster

  } else if (clustering_method == "leiden") {

    knn_graph <- RANN::nn2(clustering_input, k = graph_distance)$nn.idx
    edges <- do.call(rbind, lapply(1:nrow(knn_graph), function(i) {
      data.frame(from = i, to = knn_graph[i, ])
    }))
    g <- igraph::graph_from_data_frame(edges, directed = FALSE)

    result_clusters <- leidenbase::leiden_find_partition(g, resolution_parameter = leiden_resolution)
    CygnusObject@ev_meta[['leiden_clusters']] <- as.integer(result_clusters$membership)

  } else {
    stop("Invalid 'clustering_method'. Choose one of: 'kmeans', 'leiden', or 'hdbscan'.")
  }

  return(CygnusObject)

}

