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

#' Create Cygnus by Launching Shiny UI for Column Selection
#'
#' This function launches a Shiny app that allows users to select columns for markers
#' and metadata interactively. It then returns the selected columns.
#'
#' @param data.path Character string specifying the path to the data file.
#' @return A list containing selected columns for markers and metadata.
#' @export
specifyColumns <- function(data.path) {
  result <- shiny::shinyApp(

    # user interface
    ui = shiny::fluidPage(
      shiny::titlePanel("Select Columns for Cygnus Object"),
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::uiOutput("column_selectors"),
          shiny::actionButton("submit", "Create Cygnus Object")
        ),
        shiny::mainPanel(
          shiny::verbatimTextOutput("summary")
        )
      )
    ),

    # server
    server = function(input, output, session) {
      data <- shiny::reactiveVal(NULL)

      shiny::observe({
        req(data.path)
        data(read.csv(data.path))
      })

      output$column_selectors <- shiny::renderUI({
        req(data())
        df <- data()
        list(
          shiny::selectInput("markers_col", "Select Columns for Markers",
                             choices = names(df), selected = NULL, multiple = TRUE),
          shiny::selectInput("meta_col", "Select Columns for Metadata",
                             choices = names(df), selected = NULL, multiple = TRUE)
        )
      })

      shiny::observeEvent(input$submit, {
        req(data())
        markers_col <- input$markers_col
        meta_col <- input$meta_col

        shiny::stopApp(list(markers_col = markers_col, meta_col = meta_col))
      })
    }
  )
  selected_vals <- shiny::runApp(result)
  return(selected_vals)
}



#' Create a Cygnus Object
#'
#' This function creates a Cygnus Object from the AIVIA output file format.
#' If no columns are provided, a Shiny app will be launched to allow
#' interactive selection of columns.
#'
#' @param data.path Character string specifying the path to the data file.
#' @param markers_col A vector of column names or indices representing markers.
#' @param meta_col A vector of column names or indices representing metadata.
#' @return An object of class \code{CygnusObject}.
#' @export
CreateCygnus <- function(data.path, markers_col = NULL, meta_col = NULL) {
  if (is.null(markers_col) || is.null(meta_col)) {
    cols <- specifyColumns(data.path)
    markers_col <- cols$markers_col
    meta_col <- cols$meta_col
  }

  data <- read.csv(data.path)

  obj <- new("CygnusObject",
             matrices = list('Raw_Score' = data[markers_col]),
             ev_meta = data[meta_col],
             markers_meta = list(),
             marker_analysis = NULL,
             dim_red = list())

  return(obj)
}
