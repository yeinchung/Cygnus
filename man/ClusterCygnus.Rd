% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering.R
\name{ClusterCygnus}
\alias{ClusterCygnus}
\title{Perform Clustering on Dimensionality Reduced Data}
\usage{
ClusterCygnus(
  CygnusObject,
  use.dims = "PCA",
  use.markers = TRUE,
  relevant_markers = TRUE,
  matrix_name = "scaled_exp_matrix",
  n_clusters = 10
)
}
\arguments{
\item{CygnusObject}{A `CygnusObject` containing the data matrices and metadata for analysis.}

\item{use.dims}{A character string specifying the dimensionality reduction technique to use.
Valid options are "PCA" (default), "tSNE", or "UMAP".}

\item{use.markers}{A logical value indicating whether to use marker genes for clustering (default is TRUE).}

\item{relevant_markers}{A logical value indicating whether to use only relevant markers for clustering
(default is TRUE).}

\item{matrix_name}{A character string specifying the name of the matrix to use from the `CygnusObject`.
The default is "scaled_exp_matrix".}

\item{n_clusters}{An integer specifying the number of clusters to form (default is 10).}
}
\value{
A `CygnusObject` with the cluster assignments added to the `ev_meta` slot.
}
\description{
Clusters data in the `CygnusObject` using specified dimensionality reduction techniques and K-means clustering.
The function updates the `CygnusObject` with the resulting cluster assignments.
}
\details{
The function performs the following steps:
1. Extracts the specified expression matrix from the `CygnusObject`.
2. Filters the matrix to include only relevant markers if specified.
3. Retrieves the selected dimensionality reduction results (`PCA`, `tSNE`, or `UMAP`).
4. Uses the reduced dimensions for K-means clustering to assign each sample to a cluster.
5. Stores the cluster assignments in the `ev_meta` slot of the `CygnusObject`.
}
