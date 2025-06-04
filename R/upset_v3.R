# Load necessary libraries
library(ggplot2)      # For data visualization
library(dplyr)        # For data manipulation
library(tidyverse)    # For data wrangling and visualization
library(ComplexUpset) # For creating complex upset plots
library(Cygnus)       # For single-cell data analysis


# Define the path to the data file
data.path <- "/Users/yeinchung/Desktop/SEAM_HumanPatient_20240611 (1).csv"
all_data <- read.csv(data.path)

write.csv(all_data[all_data$stage %in% c("F1", "F2", "M1"), ], "/Users/yeinchung/Desktop/cygnus/healthy.csv")
write.csv(all_data[all_data$stage %in% c("T1-1", "T1-2", "T1-3"), ], "/Users/yeinchung/Desktop/cygnus/t1.csv")
write.csv(all_data[all_data$stage %in% c("T2-1", "T2-2", "T2-3"), ], "/Users/yeinchung/Desktop/cygnus/t2.csv")
write.csv(all_data[all_data$stage %in% c("T3-1", "T3-2", "T4-1"), ], "/Users/yeinchung/Desktop/cygnus/t3_t4.csv")


relevant <- "/Users/yeinchung/Desktop/cygnus/healthy.csv"

# Define the markers of interest
#raw_markers <-  c("PanEV", "EpCAM", "MET",
#                 "SDC1", "EGFR", "ADAM10",
#                "CTSH", "PDL1", "HER2")

# Create a Cygnus object using the data file and specified markers
cyg <- CreateCygnus(relevant,
                    meta_col = "stage", markers_col = c("EpCAM", "MET", "SDC1", "EGFR", "SP.B", "CTSH", "PDL1", "ADAM10", "MUC1","HER2"))



# Extract a binary expression matrix based on a threshold of 1000
cyg <- Cygnus::createBinaryMatrix(cyg, thresholds = 0.5)
binary_exp_matrix <- cyg@matrices[["binary_exp_matrix"]][, c("EpCAM", "MET", "SDC1", "EGFR", "SP.B", "CTSH", "PDL1", "ADAM10", "MUC1","HER2")]


# Prepare a membership matrix by combining binary expression data
membership_matrix <- as.data.frame(binary_exp_matrix) %>%
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
  count(Intersection, name = "Observed")


# Simulate null distribution of deviations by shuffling marker expressions
iteration_no = 100 # Number of simulations (1000 for final analysis)

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
    count(Intersection, name = "Shuffled_Observed")

  # Calculate deviations between shuffled and expected counts
  deviation <- expected_prob_df %>%
    left_join(shuffled_counts, by = "Intersection") %>%
    mutate(
      Shuffled_Observed = replace_na(Shuffled_Observed, 0), # Replace NA with 0
      Deviation = Shuffled_Observed - Expected_Count   # Compute deviation
    ) %>%
    select(Intersection, Deviation)

  deviation
}, simplify = FALSE)


# Aggregate null deviations across simulations
null_deviation_summary <- bind_rows(null_deviation_distributions, .id = "Simulation") %>%
  group_by(Intersection) %>%
  summarise(
    MeanDeviation = mean(Deviation, na.rm = TRUE),  # Mean deviation
    SDDeviation = sd(Deviation, na.rm = TRUE),      # Standard deviation of deviation
    .groups = 'drop'
  )

# Calculate deviation statistics for the actual observed data
deviation_stats <- observed %>%
  left_join(expected_prob_df, by = "Intersection") %>%
  mutate(
    Actual_Deviation = Observed - Expected_Count     # Compute actual deviation
  ) %>%
  left_join(null_deviation_summary, by = "Intersection") %>%
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
  left_join(deviation_stats, by = "Intersection")


# Prepare data for upset plot
membership_matrix_mixed <- membership_matrix_long %>%
  mutate(
    significant = ifelse(is.na(p_adj), "",  # Add significance labels
                         ifelse(p_adj < 0.005, "*", ""))
  )



# intersection count, min_degree 설정
threshold_count<-1
min_degree_setting <-2

#  Intersection 별 count 계산 (treshold_count 이상만 유지하기 위함 )
intersection_counts <- membership_matrix_mixed %>%
  count(Intersection, name = "Observed") %>%  #  각 cell 갯수  계산
  filter(Observed >= threshold_count) #  1threshold_count 미만 조합 제외

# 유효한 intersection 만 따로 저장
valid_intersections <- intersection_counts$Intersection

# 전체 data에서 threshold_count 이상 조합만 남긴 새로운 데이터
filtered_matrix <- membership_matrix_mixed %>%
  filter(Intersection %in% valid_intersections)


# 💡 각 row에 대해 몇 개의 마커가 양성인지(degree) 계산
marker_matrix <- filtered_matrix[, colnames(binary_exp_matrix)]  #  마커 값만 추출
filtered_matrix$Degree <- rowSums(marker_matrix)  #  각 셀마다 양성 마커 수 합산


# min_degree 이상인 데이터만 선택  (upset 그릴 때 실제 표시 대상이 됨)
plot_data <- filtered_matrix %>%
  filter(Degree >= min_degree_setting)

#  위에서 필터링한 조합 중 가장 많은 intersection count를 구함 (y축 최대값 계산용)
max_count <- plot_data %>%
  count(Intersection) %>%
  pull(n) %>%
  max()

# upset plot
upset_plot <- upset(
  plot_data,
  intersect = colnames(binary_exp_matrix),
  min_degree = min_degree_setting,
  name = "biomarkers",
  stripes = c('white', 'grey90'),  # Color stripes for the plot #deepskyblue1
  set_sizes = (
    upset_set_size()       # Add set sizes to the plot
    + geom_text(
      aes(label = ..count..),  # Add counts to set sizes
      stat = "count",
      hjust = -0.5,
      size = 3,
      color = "white"
    )
    + theme(
      axis.text.x = element_blank(),   # Remove x-axis text
      axis.ticks.x = element_blank(),  # Remove x-axis ticks
    )
  ),


  base_annotations=list(
    'Intersection size'=intersection_size(
      counts= TRUE,         # Show intersection sizes
      bar_number_threshold = 1,   # Show all bars
      mapping =aes(fill = Actual_Deviation),   # Fill bars by actual deviation
      text_colors=c(
        on_background='brown', on_bar='yellow' # Text colors
      ),
    )
    +
      annotate(   # Add annotation for significance
        geom='text', x=Inf, y=Inf,
        label = paste('* : p-value < 0.005' ),
        vjust=1.5, hjust=1,
      )

    +
      geom_text(  # Add significance stars to bars
        aes(
          label = significant,
          y = Observed + (Observed*0.05)
        ),
        size = 5,
        vjust = -0.5,
        color = "black"
      )

    +
      scale_fill_gradient2(  # Color gradient for deviations
        low = "orange", mid = "pink", high = "darkblue",
        midpoint = 0, na.value = "gray",
        guide = guide_colorbar(title = "Deviation")
      )

    + ylab('Intersection size')  # Y-axis label
    + scale_y_continuous(limits = c(0, 1500))    ), # it was limits = c(0, max_count * 1.2) for automatic scaling

  width_ratio = 0.1
)

print(upset_plot)



