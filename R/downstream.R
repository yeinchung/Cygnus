#' Visualize Proportions in Bar Plot
#'
#' This function creates a bar plot to visualize the proportions of a specified
#' group within the data. The bars are colored according to a specified metadata
#' column and can be grouped by another metadata column.
#'
#' @param data A Cygnus object containing EV metadata and matrices.
#' @param color_by A column name from the EV metadata to use for coloring the proportions in the plot. This defines the different groups/colors in the bar plot.
#' @param split_by A column name from the EV metadata that defines how to group the data on the x-axis.
#'
#' @return A ggplot2 object representing the bar plot of proportions.
#' @export
plotBar <- function(data, color_by = NULL, split_by = NULL){
  pt <-table(data@ev_meta[[color_by]], data@ev_meta[[split_by]]) # include error
  pt <- as.data.frame(pt)

  color_palette <- RColorBrewer::brewer.pal(length(levels(data@ev_meta[[split_by]])), "Set2")

  ggplot2::ggplot(pt, ggplot2::aes(x = Var2, y = Freq, fill = Var1)) +
    ggplot2::theme_bw(base_size = 15) +
    ggplot2::geom_col(position = "fill", width = 0.5) +
    ggplot2::xlab(split_by) +
    ggplot2::ylab("Proportion") +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::scale_fill_manual(values = color_palette)
}



#' Visualize Average Expression as a Line Graph
#'
#' This function generates a line plot to visualize the average expression of specified markers
#' across different groups. The x-axis represents the groups defined by a metadata column,
#' and the y-axis shows the mean expression values for each marker.
#'
#' @param data A Cygnus object containing EV metadata and matrices.
#' @param x_axis A column name from the EV metadata to use for the x-axis. This column should define the groups for the line plot.
#' @param plot_markers A vector of marker names to plot. If set to "ALL", all markers will be included in the plot. Default is "ALL".
#' @param fold_change Logical; if TRUE, the fold change will be calculated for the markers. Default is TRUE.
#' @param matrix_name The name of the matrix within the Cygnus object to use for the plot. Default is "Raw_Score".
#' @param order_x_axis A vector specifying the order of levels for the x-axis. If NULL, the default order will be used. Default is NULL.
#'
#' @return A ggplot2 object representing the line plot of average expression values.
#' @importFrom dplyr %>%
#' @export
plotLine <- function(data, x_axis, plot_markers = "ALL", fold_change=T, matrix_name="Raw_Score", order_x_axis=NULL){
  m <- data@matrices[[matrix_name]]
  if(plot_markers != "ALL"){m <- m[, plot_markers]}
  group <- data@ev_meta[[x_axis]]
  m <- cbind(m, group)
  if(is.null(order_x_axis)){
    m$group <- factor(m$group)
  }else{
      m$group <- factor(m$group, levels = order_x_axis)
  }
  m.small <- m %>% group_by(group) %>% summarize(across(where(is.numeric), mean, na.rm = TRUE))

  m_long <- m.small %>%
    pivot_longer(cols = -group, names_to = "variable", values_to = "mean_val")

  ggplot(data = m_long, aes(x = group, y = mean_val, colour = variable)) +
    geom_line(aes(group = variable)) +
    geom_point() +
    labs(title = paste0("Average intensity scores from ", matrix_name),
         x = x_axis,
         y = "Expression",
         colour = "Markers") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

}
