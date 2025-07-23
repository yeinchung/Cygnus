#' Generate an UpSet Plot for Relevant Markers
#'
#' @param data An object of class \code{CygnusObject} containing the binary expression matrix.
#' @param markers Character vector specifying which markers to include in the UpSet plot. If NULL, all relevant markers are used.
#' @param ... Additional arguments to pass to the upset function.
#' @return An UpSet plot visualizing the intersections of relevant markers.
plotUpSet_basic <- function(data, markers = NULL, ...) {
  # Check if binary_exp_matrix exists
  if (!"binary_exp_matrix" %in% names(data@matrices)) {
    stop("Binary expression matrix not found. Please run createBinaryMatrix first.")
  }

  binary_matrix <- data@matrices$binary_exp_matrix

  relevant_markers <- data@markers_meta$marker[data@markers_meta$relevant]
  if (is.null(relevant_markers) || length(relevant_markers) == 0) {
    stop("No relevant markers found.")
  }

  if (!is.null(markers)) {
    relevant_markers <- relevant_markers[relevant_markers %in% markers]
    if (length(relevant_markers) == 0) {
      stop("None of the specified markers were found in relevant markers.")
    }
  }

  binary_matrix <- binary_matrix[, relevant_markers, drop = FALSE]

  binary_df <- as.data.frame(binary_matrix)
  colnames(binary_df) <- relevant_markers

  UpSetR::upset(binary_df, ...)
}


#' Generate Marker Colocalization Deviation and Statistics
#' @param data An object of class \code{CygnusObject} containing the binary expression matrix.
#' @param matrix_name The matrix within the data object that contains the binary expression matrix. The default is set to binary_exp_mat
#' @param iteration_no The number of iterations for simulation
#' @return A dataframe of deviation and statistics of each intersection
#' @export
getColocalizedMarkers <- function(data, matrix_name = "binary_exp_matrix", iteration_no = 100){
  # Check if binary_exp_matrix exists
  if (!"binary_exp_matrix" %in% names(data@matrices)) {
    stop("Binary expression matrix not found. Please run createBinaryMatrix first.")
  }

  binary_matrix <- data@matrices$binary_exp_matrix

  relevant_markers <- data@markers_meta$marker[data@markers_meta$relevant]
  if (is.null(relevant_markers) || length(relevant_markers) == 0) {
    stop("No relevant markers found.")
  }

  binary_matrix <- binary_matrix[, relevant_markers, drop = FALSE]

  binary_df <- as.data.frame(binary_matrix)
  colnames(binary_df) <- relevant_markers

  print("Calculating expected and observed probabilities...")

  # Prepare a membership matrix by combining binary expression data
  membership_matrix <- as.data.frame(binary_df) %>%
    mutate(Intersection = apply(., 1, paste0, collapse = ""))  # Create a unique identifier for each combination of markers

  # Calculate the total number of observations (cells)
  total_observations <- nrow(membership_matrix)


  # Calculate the probability of each marker being expressed (1)
  prob <- membership_matrix %>%
    summarise(across(everything(), ~ mean(. == 1))) # Compute the mean expression for each marker
  within(prob, rm(Intersection))  # Remove the Intersection column for further calculations


  # Compute expected probabilities for each intersection under the null hypothesis (independent expression)
  present_intersections <- unique(membership_matrix$Intersection)   # Unique combinations of markers
  expected_prob_df <- data.frame(Intersection = present_intersections, Expected_Probability = NA)


  prob <- within(prob, rm(Intersection))
  # Calculate expected probabilities for each intersection
  for (i in seq_along(present_intersections)) {
    intersection <- present_intersections[i]
    temp_prob <- 1

    for (j in seq_along(prob)) {
      marker_prob <- prob[[j]]  # Probability of the j-th marker being expressed
      marker_presence <- as.integer(substr(intersection, j, j))  # Presence (1) or absence (0) of the marker in the intersection

      if (marker_presence == 1) {
        temp_prob <- temp_prob * marker_prob   # Multiply by the probability of presence
      } else {
        temp_prob <- temp_prob * (1 - marker_prob)  # Multiply by the probability of absence
      }
    }
    # Store the expected probability for the intersection
    expected_prob_df$Expected_Probability[i] <- temp_prob
  }


  # Convert expected probabilities to expected counts
  expected_prob_df <- expected_prob_df %>%
    mutate(Expected_Count = Expected_Probability * total_observations)

  # Calculate observed counts for each intersection from the actual data
  observed <- membership_matrix %>%
    dplyr::count(Intersection, name = "Observed")

  print("Running simulation to calculate the null distribution of deviations...")

  # Simulate null distribution of deviations
  set.seed(1009)           # Set seed for reproducibility
  null_deviation_distributions <- replicate(iteration_no, {
    # Shuffle the marker expressions to simulate random co-expression
    shuffled <- apply(membership_matrix[, -ncol(membership_matrix)], 2, sample)

    # Create Intersection column from shuffled data
    shuffled_long <- as.data.frame(shuffled) %>%
      mutate(Intersection = apply(., 1, paste0, collapse = ""))

    # Count intersections in shuffled data
    shuffled_counts <- shuffled_long %>%
      dplyr::count(Intersection, name = "Shuffled_Observed")

    # Calculate deviations between shuffled and expected counts
    deviation <- expected_prob_df %>%
      dplyr::left_join(shuffled_counts, by = "Intersection") %>%
      mutate(
        Shuffled_Observed = tidyr::replace_na(Shuffled_Observed, 0), # Replace NA with 0
        Deviation = Shuffled_Observed - Expected_Count   # Compute deviation
      ) %>%
      select(Intersection, Deviation)

    deviation
  }, simplify = FALSE)


  # Aggregate null deviations across simulations
  null_deviation_summary <- dplyr::bind_rows(null_deviation_distributions, .id = "Simulation") %>%
    group_by(Intersection) %>%
    summarise(
      MeanDeviation = mean(Deviation, na.rm = TRUE),  # Mean deviation
      SDDeviation = sd(Deviation, na.rm = TRUE),      # Standard deviation of deviation
      .groups = 'drop'
    )

  # Calculate deviation statistics for the actual observed data
  deviation_stats <- observed %>%
    dplyr::left_join(expected_prob_df, by = "Intersection") %>%
    mutate(
      Actual_Deviation = Observed - Expected_Count     # Compute actual deviation
    ) %>%
    dplyr::left_join(null_deviation_summary, by = "Intersection") %>%
    mutate(
      p_value = pnorm(Actual_Deviation / SDDeviation, lower.tail = FALSE),   # Compute p-value
      p_adj = p.adjust(p_value, method = "bonferroni")    # Adjust p-values for multiple testing
    )


  # Cap p-values for visualization and avoid -Inf in log transformation
  deviation_stats <- deviation_stats %>%
    mutate(
      p_adj_capped = pmax(p_adj, 1e-300),    # Cap p-values at a minimum threshold
      log10_p_adj = -log10(p_adj_capped)     # Compute -log10 adjusted p-values
    )

  # Merge deviation statistics back into the membership matrix for plotting
  membership_matrix_long <- membership_matrix %>%
    dplyr::left_join(deviation_stats, by = "Intersection")


  # Prepare data for upset plot
  membership_matrix_mixed <- membership_matrix_long %>%
    mutate(
      significant = ifelse(is.na(p_adj), "",  # Add significance labels
                           ifelse(p_adj < 0.005, "*", ""))
    )

  return(as.data.frame(membership_matrix_mixed))
  print(head(membership_matrix_mixed))
}

#' Generate an UpSet Plot from Colocalization Statistics
#'
#' @param data An object of class \code{CygnusObject} containing the binary expression matrix.
#' @param colocal_data A data frame of deviations and statistics from running getColocalizedMarkers
#' @param matrix_name The matrix within the data object that contains the binary expression matrix. The default is set to binary_exp_mat
#' @param threshold_count An integer indicating the minimum number of elements in an intersection
#' @param min_degree_setting An integer indicating the minimum degree of intersection for plotting
#' @param upper_limit Specifies maximum value for y-axis
#' @export
plotUpset <- function(data, colocal_data, matrix_name = "binary_exp_mat", threshold_count = 1, min_degree_setting = 1){
  if (!"binary_exp_matrix" %in% names(data@matrices)) {
    stop("Binary expression matrix not found. Please run createBinaryMatrix first.")
  }

  binary_matrix <- data@matrices$binary_exp_matrix

  relevant_markers <- data@markers_meta$marker[data@markers_meta$relevant]
  if (is.null(relevant_markers) || length(relevant_markers) == 0) {
    stop("No relevant markers found.")
  }

  binary_matrix <- binary_matrix[, relevant_markers, drop = FALSE]

  binary_df <- as.data.frame(binary_matrix)
  colnames(binary_df) <- relevant_markers


  #  Intersection 별 count 계산 (treshold_count 이상만 유지하기 위함 )
  intersection_counts <- colocal_data %>%
    dplyr::count(Intersection, name = "Observed") %>%  #  각 cell 갯수  계산
    filter(Observed >= threshold_count) #  1threshold_count 미만 조합 제외

  # 유효한 intersection 만 따로 저장
  valid_intersections <- intersection_counts$Intersection

  # 전체 data에서 threshold_count 이상 조합만 남긴 새로운 데이터
  filtered_matrix <- colocal_data %>%
    filter(Intersection %in% valid_intersections)


  # 💡 각 row에 대해 몇 개의 마커가 양성인지(degree) 계산
  marker_matrix <- filtered_matrix[, colnames(binary_df)]  #  마커 값만 추출
  filtered_matrix$Degree <- rowSums(marker_matrix)  #  각 셀마다 양성 마커 수 합산


  # min_degree 이상인 데이터만 선택  (upset 그릴 때 실제 표시 대상이 됨)
  plot_data <- filtered_matrix %>%
    filter(Degree >= min_degree_setting)

  #  위에서 필터링한 조합 중 가장 많은 intersection count를 구함 (y축 최대값 계산용)
  max_count <- plot_data %>%
    dplyr::count(Intersection) %>%
    dplyr::pull(n) %>%
    max()

  # Get the maximum set size.
  # This value is to be uese to scale the set size plot
  max_set_size <- plot_data[, relevant_markers] %>%
    colSums() %>%
    max()

  # upset plot
  upset_plot <- ComplexUpset::upset(
    plot_data,
    intersect = colnames(binary_df),
    min_degree = min_degree_setting,
    name = "Biomarker combinations",
    stripes = c('white'),  # Color stripes for the plot #deepskyblue1
    #stripes = c('white', 'grey90'),  # Color stripes for the plot #deepskyblue1
    set_sizes = (
      ComplexUpset::upset_set_size(
        geom=geom_bar(width=0.8, stat='count'),
        position='right')
      + geom_text(
        #aes(label = ..count..),  # Add counts to set sizes
        aes(label = after_stat(count)),
        stat = "count",
        hjust = -0.2,
        vjust = 0.5,
        size.unit = "pt",
        size = 9,
        #color = "white"
        #color = 'blue'
        color = "black"
      )
      + ylab('Set size')
      + scale_y_continuous(limits = c(0, max_set_size*2),
                           expand = expansion(mult = c(0, 0.1))
      )
      + theme(
        axis.title.x = element_text(size = 11),
        axis.text.x = element_blank(),   # Remove x-axis text
        axis.ticks.x = element_blank(),  # Remove x-axis ticks
      )
    ),
    guides='over',

    base_annotations=list(
      'Intersection size'=ComplexUpset::intersection_size(
        counts= TRUE,         # Show intersection sizes
        bar_number_threshold = 1,   # Show all bars
        mapping =aes(fill = Actual_Deviation),   # Fill bars by actual deviation
        text_colors=c(
          on_background='brown', on_bar='yellow' # Text colors
        ),
        text <- list(size.unit = "pt", size = 9),
      )
      +
        annotate(   # Add annotation for significance
          geom='text', x=Inf, y=Inf,
          label = substitute(paste("*", italic('p'), " < 0.005")),
          size.unit = "pt", size = 10,
          #vjust=1.5, hjust= 1.,
          vjust=3, hjust= 1.5)

      +
        geom_text(  # Add significance stars to bars
          aes(
            #label = significant,
            label = ifelse(p_adj < 0.005, "*", ""),
            #y = Observed + (Observed*0.05)
            y = Observed  # Offset text above the bar
          ),
          vjust = -1, hjust = 0.5,
          size.unit = "pt", size = 11,
          color = "black",
          #color = "red"
        )
      +
        scale_fill_gradient2(  # Color gradient for deviations
          low = "orange", mid = "pink", high = "darkblue",
          midpoint = 0, na.value = "gray",
          # make guide_colorbar title font size larger
          guide = guide_colorbar(
            title = "Deviation",
            title.theme = element_text(size = 10),
            theme = theme(legend.text = element_text(size = 9),
                          legend.key.width  = unit(0.5, "lines"),
                          legend.key.height = unit(5, "lines"))
          )
        )

      + ylab('Intersection size')  # Y-axis label
      + theme(axis.title.y = element_text(size = 12, vjust = -6),
              axis.text.y = element_text(size = 11))
      #      + scale_y_continuous(limits = c(0, 1500))     # it was limits = c(0, max_count * 1.2) for automatic scaling
      + scale_y_continuous(limits = c(0, max_count * 1.2), expand = expansion(mult = c(0, 0.1))) # Add some space at the top
      + theme(axis.line = element_line(color = "black", size = 0.25),  # Set axis line color and size
              axis.ticks = element_line(color = "black", size = 0.25),  # Set axis ticks color and size
              panel.grid.major = element_blank(),  # Remove major grid lines
              panel.grid.minor = element_blank(),  # Remove minor grid lines
              panel.background = element_rect(fill = "#FFFEE6", color = "white"),  # Set panel background color
              #panel.background = element_blank(),  # Remove background color
              #plot.margin = margin(10, 10, 10, 10)  # Set plot margins
      )
    ),
    width_ratio = 0.15,
    height_ratio = 0.8,
    themes = ComplexUpset::upset_default_themes(
      #text = element_text(color= 'red'),
      axis.text = element_text(size = 11),  # Set axis text size
    )  # Set text size for the plot
  )
  print(upset_plot)  # Print the plot with specified dimensions
  return(upset_plot)

}

#' Generate an Volcano Plot from Colocalization Statistics
#'
#' @param colocal_data A data frame of deviations and statistics from running getColocalizedMarkers
#' @param sig_threshold The threshold for significance
#' @export
plotVolcano <- function(colocal_data, sig_threshold = 0.05){
  significance_threshold <- 0.05
  log_threshold <- -log10(significance_threshold)

  deviation_stats <- unique(colocal_data)

  ggplot(deviation_stats, aes(x = Actual_Deviation, y = log10_p_adj, label = Intersection)) +
    geom_point(aes(color = log10_p_adj, size = abs(Actual_Deviation)), alpha = 0.7) +
    ggrepel::geom_text_repel(data = subset(subset(deviation_stats, log10_p_adj > 30), Actual_Deviation > 10),
                    aes(label = Intersection),
                    size = 4, fontface = "bold", color = "darkblue",
                    box.padding = 0.6, max.overlaps = 15) +
    geom_hline(yintercept = log_threshold, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_gradientn(colors = c("blue", "purple", "red"),
                          limits = c(0, max(deviation_stats$log10_p_adj))) +
    labs(
      x = "Actual Deviation",
      y = expression(-log[10](adjusted~p~value)),
      color = "-log10(p_adj)",
      size = "Deviation Magnitude"
    ) +
    theme_classic(base_size = 15) +
    theme(legend.position = "right")
}
