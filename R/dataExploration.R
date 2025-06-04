#' Plot Distributions of Selected Markers
#'
#' This function plots the distribution of expression levels for selected markers in an expression matrix.
#' The user can specify a subset of markers to plot, or plot all markers by default.
#' Additionally, the user can group histograms by a specified column in ev_meta.
#'
#' @param obj An object of class \code{CygnusObject} containing the expression matrix and metadata.
#' @param plot_markers A vector of marker names to plot. If set to "ALL", distributions for all markers will be plotted. Default is "ALL".
#' @param matrix Character string specifying the name of the matrix to use. Default is "Raw_Score".
#' @param group_by A character string specifying a column in \code{ev_meta} to group histograms by. Default is \code{NULL} (no grouping).
#' @param group_colors A named vector of colors for groups. Names should correspond to unique values in the group_by column. If not specified, defaults to the RColorBrewer Set2 palette.
#' @return A series of histogram plots, each representing the distribution of expression levels for the specified markers.
#' @export
plotDistribution <- function(
    obj,
    plot_markers = "ALL",
    matrix = "Raw_Score",
    group_by = NULL,
    group_colors = NULL
) {
  # Check if the specified matrix exists
  if (!(matrix %in% names(obj@matrices))) {
    stop(paste("Matrix", matrix, "not found in CygnusObject"))
  }

  matrix_data <- obj@matrices[[matrix]]

  # Check if plot_markers is "ALL" or a subset
  if (plot_markers == "ALL") {
    plot_markers <- colnames(matrix_data)
  } else {
    plot_markers <- intersect(plot_markers, colnames(matrix_data))
  }

  if (length(plot_markers) == 0) {
    stop("No markers to plot.")
  }

  # Check if group_by is provided and exists in ev_meta
  if (!is.null(group_by)) {
    if (!(group_by %in% names(obj@ev_meta))) {
      stop(paste("Group column", group_by, "not found in ev_meta"))
    }
    # Extract unique groups, excluding NAs
    groups <- unique(na.omit(obj@ev_meta[[group_by]]))

    # Ensure groups are valid
    if (length(groups) == 0) {
      stop("No valid groups found in the group_by column.")
    }

    # Assign colors to groups
    if (is.null(group_colors)) {
      brewer_colors <- RColorBrewer::brewer.pal(max(3, length(groups)), "Set2")
      group_colors <- setNames(brewer_colors[1:length(groups)], groups)
    } else {
      # Check if provided colors match the groups
      if (!all(groups %in% names(group_colors))) {
        stop("Not all groups have corresponding colors in group_colors.")
      }
    }
  } else {
    groups <- NULL
  }

  num_cols <- length(plot_markers)
  num_rows <- max(1, ceiling(num_cols / 3))

  graphics::par(mfrow = c(num_rows, min(3, num_cols)))
  graphics::par(mar = c(2, 2, 2, 1))

  for (marker in plot_markers) {
    if (!is.null(groups)) {
      # Plot separate histograms for each group with different colors
      for (group in groups) {
        group_indices <- obj@ev_meta[[group_by]] == group
        hist(matrix_data[group_indices, marker], breaks = 200,
             main = paste(marker, "-", group),
             xlab = "", ylab = "",
             col = group_colors[group], border = group_colors[group])
      }
    } else {
      # Plot a single histogram
      hist(matrix_data[, marker], breaks = 1000,
           main = paste(marker),
           xlab = "", ylab = "")
    }
  }

  graphics::par(mfrow = c(1, 1))
}



#' Plot Average Expression Heatmap by Group
#'
#' This function generates a heatmap of average marker expressions, grouped by a specified metadata column.
#' The heatmap allows customization of clustering distance, color palette, font size, and scaling method.
#'
#' @param data An object of class \code{CygnusObject} containing an expression matrix and metadata.
#' @param group_column Character string specifying the metadata column to group by.
#' @param clustering_distance Character string specifying the distance metric for clustering columns. Default is "euclidean".
#' @param colors A color palette to use for the heatmap. Default is a reversed "RdYlBu" palette from RColorBrewer.
#' @param fontsize Numeric value specifying the font size for heatmap text. Default is 8.
#' @param scale Character string specifying whether the data should be scaled by 'row', 'column', or 'none'. Default is 'column'.
#' @param cluster_rows Logical value indicating whether to cluster rows in the heatmap. Default is FALSE.
#' @param na.rm Logical value indicating whether to remove NA values before calculating averages. Default is TRUE.
#' @return A heatmap plot showing average marker expressions grouped by the specified metadata column.
#'
#' @importFrom dplyr %>%
#'
#' @export
plotAvgHeatmap <- function(data, group_column,
                           clustering_distance = "euclidean",
                           colors = rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = 11, "RdYlBu"))(100)),
                           fontsize = 8,
                           scale = 'column',
                           cluster_rows = FALSE,
                           na.rm = TRUE) {

  expression_matrix <- data@matrices$Raw_Score
  group_metadata <- data@ev_meta[[group_column]]

  combined_data <- cbind(expression_matrix, group_metadata)
  avg_marker_expressions <- combined_data %>%
    dplyr::group_by(group_metadata) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), mean, na.rm = na.rm)) %>%
    dplyr::select(-group_metadata) %>%
    as.matrix()

  pheatmap::pheatmap(avg_marker_expressions,
                    # clustering_distance_cols = clustering_distance,
                     cluster_rows = cluster_rows,
                     color = colors,
                    # main = paste("Average Marker Expressions by", group_column),
                     fontsize = fontsize,
                     scale = scale,
                     labels_row = unique(group_metadata))

}

#' Scale Expression Matrix by Maximum Value
#'
#' This function scales the expression matrix by dividing each marker's values by the maximum value in that marker.
#' The scaled matrix is added as an additional layer in the data object.
#'
#' @param data An object of class \code{CygnusObject} containing an expression matrix.
#' @param matrix_name Character string specifying the name for the new scaled matrix layer. Default is "scaled_exp_matrix".
#' @return The data object with an additional layer containing the scaled expression matrix.
#' @export
scaleExpressionMatrix <- function(data, matrix_name = "scaled_exp_matrix") {

  expression_matrix <- data@matrices$Raw_Score

  max_values <- apply(expression_matrix, 2, max, na.rm = TRUE)
  scaled_matrix <- sweep(expression_matrix, 2, max_values, `/`)

  data@matrices[[matrix_name]] <- scaled_matrix

  return(data)
}

#' Create Binary Expression Matrix Based on Thresholds
#'
#' This function creates a binary matrix where each entry is 1 if the marker's raw score exceeds a specified threshold,
#' and 0 otherwise. The thresholds can be customized for each marker using a vector. If no vector is provided,
#' a default threshold of 100 is used for all markers. The binary matrix is added as an additional layer in the data object.
#'
#' @param data An object of class \code{CygnusObject} containing an expression matrix.
#' @param thresholds Numeric vector specifying the thresholds for each marker. If a single number is provided, it is used for all markers. Default is 100.
#' @param matrix_name Character string specifying the name for the new binary matrix layer. Default is "binary_exp_matrix".
#' @return The data object with an additional layer containing the binary expression matrix.
#' @export
createBinaryMatrix <- function(data, thresholds = 100, matrix_name = "binary_exp_matrix") {

  if (!("Raw_Score" %in% names(data@matrices))) {
    stop("Raw_Score matrix not found in CygnusObject")
  }

  expression_matrix <- data@matrices$Raw_Score
  num_markers <- ncol(expression_matrix)

  if (length(thresholds) == 1) {
    thresholds <- rep(thresholds, num_markers)
  } else if (length(thresholds) != num_markers) {
    stop("The length of the thresholds vector must match the number of markers in the expression matrix.")
  }

  binary_matrix <- expression_matrix
  for (i in seq_len(num_markers)) {
    binary_matrix[, i] <- ifelse(expression_matrix[, i] > thresholds[i], 1, 0)
  }

  data@matrices[[matrix_name]] <- binary_matrix

  return(data)
}

#' Mark Relevant Markers
#'
#' This function updates the markers_meta slot in a CygnusObject to mark specified markers as relevant.
#' The markers_meta must be a data frame with a 'marker' column. The function adds a 'relevant' list.
#'
#' @param data An object of class \code{CygnusObject}.
#' @param relevant_markers A character vector of marker names to be marked as relevant. If NULL, no markers are marked as relevant.
#' @return The updated \code{CygnusObject} with relevant markers marked.
#' @export
markRelevantMarkers <- function(data, relevant_markers = NULL) {
  data@markers_meta[['marker']] <- colnames(data@matrices$Raw_Score)

  names_vector <- colnames(data@matrices$Raw_Score)
  data@markers_meta[['relevant']] <- stats::setNames(rep(FALSE, length(names_vector)), names_vector)

  for(marker in colnames(data@matrices$Raw_Score)) {
    if(marker %in% relevant_markers) {
      data@markers_meta[['relevant']][marker] <- TRUE
    }
  }

  return(data)
}

#' Mark Relevant Markers Using User Interface
#'
#' This function updates the markers_meta slot in a CygnusObject to mark specified markers as relevant.
#' The markers_meta must be a data frame with a 'marker' column. The function adds a 'relevant' list.
#'
#' @param data An object of class \code{CygnusObject}.
#' @return The updated \code{CygnusObject} with relevant markers marked.
#' @export
markRelevantMarkersUI <- function(data) {
  result <- shiny::shinyApp(

    # user interface
    ui = shiny::fluidPage(
      shiny::titlePanel("Select Relevant Markers"),
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::selectInput("relevant_markers", "Select Relevant Markers",
                             choices = colnames(data@matrices$Raw_Score),
                             selected = NULL, multiple = TRUE),
          shiny::actionButton("submit", "Update Cygnus Object")
        ),
        shiny::mainPanel(
          shiny::verbatimTextOutput("summary")
        )
      )
    ),

    # server
    server = function(input, output, session) {
      markers_meta <- shiny::reactiveVal(data@markers_meta)

      shiny::observeEvent(input$submit, {
        req(input$relevant_markers)

        updated_meta <- markers_meta()
        names_vector <- colnames(data@matrices$Raw_Score)
        updated_meta[['relevant']] <- stats::setNames(rep(FALSE, length(names_vector)), names_vector)

        for (marker in input$relevant_markers) {
          updated_meta[['relevant']][marker] <- TRUE
        }

        data@markers_meta <- updated_meta

        output$summary <- shiny::renderPrint({
          list(Relevant_Markers = input$relevant_markers)
        })

        shiny::stopApp(data)
      })
    }
  )
  updated_data <- shiny::runApp(result)
  return(updated_data)
}


#' Plot Distributions of Selected Markers For Threshold Selection
#'
#' This function plots the distribution of expression levels for selected markers in an expression matrix.
#' Users can interactively select regions above the x-axis to set positive thresholds for each marker.
#'
#' @param obj An object of class \code{CygnusObject} containing the expression matrix.
#' @param plot_markers A vector of marker names to plot. If set to "ALL", distributions for all markers will be plotted. Default is "ALL".
#' @param matrix Character string specifying the name of the matrix to use. Default is "Raw_Score".
#' @return A list of selected thresholds for each marker.
#' @export
selectThreshold <- function(
    obj,
    plot_markers = "ALL",
    matrix = "Raw_Score"
) {
  # Check if the specified matrix exists
  if (!(matrix %in% names(obj@matrices))) {
    stop(paste("Matrix", matrix, "not found in CygnusObject"))
  }

  matrix_data <- obj@matrices[[matrix]]

  if (plot_markers == "ALL") {
    plot_markers <- colnames(matrix_data)
  } else {
    plot_markers <- intersect(plot_markers, colnames(matrix_data))
  }

  num_cols <- length(plot_markers)

  # Ensure valid plotting parameters
  if (num_cols == 0) {
    stop("No markers to plot.")
  }

  thresholds <- shiny::shinyApp(
    ui = shiny::fluidPage(
      shiny::titlePanel("Select Thresholds for Markers"),
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::actionButton("submit", "Submit Thresholds")
        ),
        shiny::mainPanel(
          shiny::uiOutput("plots")
        )
      )
    ),
    server = function(input, output, session) {
      selected_thresholds <- shiny::reactiveVal(list())

      output$plots <- shiny::renderUI({
        plot_outputs <- lapply(plot_markers, function(marker) {
          plotly::plotlyOutput(paste0("plot_", marker))
        })
        do.call(shiny::tagList, plot_outputs)
      })

      lapply(plot_markers, function(marker) {
        output[[paste0("plot_", marker)]] <- plotly::renderPlotly({
          p <- ggplot2::ggplot(data.frame(x = matrix_data[, marker]), ggplot2::aes(x = x)) +
            ggplot2::geom_histogram(bins = 100, fill = "blue", alpha = 0.6) +
            ggplot2::labs(title = marker, x = "", y = "Frequency")

          plotly::ggplotly(p) %>%
            plotly::layout(dragmode = "select") %>%
            plotly::event_register("plotly_selected")
        })

        shiny::observeEvent(plotly::event_data("plotly_selected", source = paste0("plot_", marker)), {
          selected_data <- plotly::event_data("plotly_selected", source = paste0("plot_", marker))
          if (!is.null(selected_data)) {
            selected_range <- range(selected_data$x)
            current_thresholds <- selected_thresholds()
            current_thresholds[[marker]] <- selected_range[2]
            selected_thresholds(current_thresholds)
          }
        })
      })

      shiny::observeEvent(input$submit, {
        shiny::stopApp(selected_thresholds())
      })
    }
  )

  selected_thresholds <- shiny::runApp(thresholds)

  return(selected_thresholds)
}

#' Normalize Expression Data by Pan-EV Marker
#'
#' This method normalizes the expression matrix in the CygnusObject by dividing
#' the expression values of all markers by the expression values of a specified
#' pan-EV marker.
#'
#' @param data An object of class \code{CygnusObject}.
#' @param pan_ev_marker Character string specifying the name of the pan-EV marker.
#' @return The updated CygnusObject with normalized expression matrix.
#' @export
normalizeByPanEV <- function(data, pan_ev_marker) {
  if (!(pan_ev_marker %in% colnames(data@matrices$Raw_Score))) {
    stop(paste("Pan-EV marker", pan_ev_marker, "not found in the expression matrix."))
  }

  expression_matrix <- data@matrices$Raw_Score
  pan_ev_values <- expression_matrix[[pan_ev_marker]]

  if (any(pan_ev_values == 0)) {
    warning("Pan-EV marker contains zero values. Normalization may lead to division by zero.")
  }

  # Initialize normalized_matrix as a data frame with the same dimensions and column names
  normalized_matrix <- data.frame(matrix(NA, nrow = nrow(expression_matrix), ncol = ncol(expression_matrix)))
  colnames(normalized_matrix) <- colnames(expression_matrix)
  rownames(normalized_matrix) <- rownames(expression_matrix)

  # Perform row-wise normalization
  for (i in 1:nrow(expression_matrix)) {
    normalized_matrix[i, ] <- as.numeric(expression_matrix[i, ]) / pan_ev_values[i]
  }

  # Assign the normalized matrix to data@matrices
  data@matrices[['normalized_exp_mat']] <- normalized_matrix

  return(data)
}



#' Add IgG distribution
#'
#' @param data An object of class \code{CygnusObject}
#' @param IgG A vector for IgG distribution
#' @param name A string specifying the name of IgG distribution to be saved in \code{CygnusObject}
#'
#' @return The updated CygnusObject with IgG distribution added
#' @export
addIgG <- function(data, IgG, name){
  data@marker_analysis <- list()
  data@marker_analysis[['name']] <- IgG

  hist(IgG, main=paste0(name, "distribution"))

  return(data)
}

