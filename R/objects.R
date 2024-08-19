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
#' @export
CygnusObject <- setClass(
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
#' This function creates a Cygnus Object from the AIVIA output file format.
#' The user must specify the columns for markers and metadata.
#'
#' @param data.path Character string specifying the path to the data file.
#' @param markers_col A vector of column names or indices representing markers.
#' @param meta_col A vector of column names or indices representing metadata.
#' @return An object of class \code{CygnusObject}.
#' @export
CreateCygnus <- function(
    data.path,
    markers_col,
    meta_col
){
  data <- read.csv(data.path)
  obj <- new("CygnusObject",
             matrices =  list('Raw_Score' = data[markers_col]),
             ev_meta = data[meta_col],
             markers_meta = list(),
             marker_analysis = NULL,
             dim_red = list())
  return(obj)
}

setMethod(
  f = "show",
  signature = "CygnusObject",
  definition = function(object) {
    cat("Summary of CygnusObject:\n")
    cat("Slot: ")
    print(names(object@matrices))
    cat("Numer of EVs: ")
    print(nrow(object@matrices$Raw_Score))
  }
)


