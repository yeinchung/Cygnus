#' Visualize Proportions
#'
#'
plotBar <- function(data, color_by = NULL, split_by = NULL){
  pt <-table(data@ev_meta[[color_by]], data@ev_meta[[split_by]]) # include error
  pt <- as.data.frame(pt)
  pt$Var1 <- as.character(pt$Var1)

  ggplot2::ggplot(pt, ggplot2::aes(x = Var2, y = Freq, fill = Var1)) +
    ggplot2::theme_bw(base_size = 15) +
    ggplot2::geom_col(position = "fill", width = 0.5) +
    ggplot2::xlab("Sample") +
    ggplot2::ylab("Proportion") +
    ggplot2::theme(legend.title = ggplot2::element_blank())
}



#' Visualize Linegraph
#'
#' Makes sense for temporal data!
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
