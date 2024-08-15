#' The Cygnus Object Class
#'
#' The CygnusObject class is a data storage class that stores expression data
#' and other related information needed for the standard Cygnus workflow.
#'
#' @slot matrices List of matrices
#' @slot ev_meta Metadata for each sample
#' @slot markers_meta Metadata for each marker
#' @slot marker_analysis Results from marker characterization analyses
#' @slot dim_red Coordinates from dimensionality reduction analysis
#'
CygnuObject <- setClass(
  Class = "CygnusObject",
  slots = list(
    matrices = "list",
    ev_meta = "list",
    markers_meta = "list",
    marker_analysis = "ANY",
    dim_red = "list"
  )
)


#' Create a Cygnus Object
#'
#' Create a Cygnus Object from the AIVIA output file format. Columns for markers
#' and metadata need to be specified.
#'
CreateCygnus <- function(
    data.path,
    markers_col,
    meta_col
){
  data <- read.csv(data.path)
  obj <- new("CygnusObject",
             matrices[['Raw_scores']] = data[markers_col],
             ev_meta = data[markers_col])
  return(obj)
}
