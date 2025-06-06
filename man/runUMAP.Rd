% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dimRed.R
\name{runUMAP}
\alias{runUMAP}
\title{Perform UMAP on a Specified Matrix and Store in Dimensionality Reduction Slot}
\usage{
runUMAP(
  data,
  matrix_name = "Raw_Score",
  use_relevant = TRUE,
  dims = 3,
  use_pcs = FALSE,
  n_pcs = 5,
  n_neighbors = 30,
  min_dist = 0.3,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{data}{An object of class \code{CygnusObject} containing an expression matrix.}

\item{matrix_name}{Character string specifying the name of the matrix to use. Default is "Raw_Score".}

\item{use_relevant}{Logical indicating whether to use only relevant markers. Default is TRUE.}

\item{dims}{Numeric value specifying the number of dimensions for UMAP. Default is 2.}

\item{use_pcs}{Logical indicating whether to use PCA results as input. Default is FALSE.}

\item{n_pcs}{Number of principal components. Used only when use_pcs = TRUE. Default is 5.}

\item{n_neighbors}{Numeric value specifying the number of neighbors for UMAP. Default is 30.}

\item{min_dist}{Numeric value specifying the minimum distance for UMAP. Default is 0.3.}

\item{verbose}{Logical value indicating whether to display a progress bar. Default is TRUE.}

\item{...}{Additional arguments to pass to the \code{uwot::umap} function.}
}
\value{
The updated CygnusObject with UMAP coordinates stored in the 'dim_red' slot.
}
\description{
Perform UMAP on a Specified Matrix and Store in Dimensionality Reduction Slot
}
