# Libraries
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ComplexUpset)

# Extract binary expression matrix
cyg <- Cygnus::createBinaryMatrix(cyg, thresholds = 1000)
binary_exp_matrix <- cyg@matrices[["binary_exp_matrix"]][, c("EpCAM", "EGFR", "SDC1", "MET", "ADAM10")]

# Prepare membership matrix
membership_matrix <- as.data.frame(binary_exp_matrix) %>%
  mutate(Intersection = apply(., 1, paste0, collapse = ""))

# Total number of observations (needed for expected counts)
total_observations <- nrow(membership_matrix)

# Calculate the probability of each marker being 1
prob <- membership_matrix %>%
  summarise(across(everything(), ~ mean(. == 1)))
within(prob, rm(Intersection))

# Compute fixed expected probabilities for each intersection
present_intersections <- unique(membership_matrix$Intersection)
expected_prob_df <- data.frame(Intersection = present_intersections, Expected_Probability = NA)

prob <- within(prob, rm(Intersection))
for (i in seq_along(present_intersections)) {
  intersection <- present_intersections[i]
  temp_prob <- 1

  for (j in seq_along(prob)) {
    marker_prob <- prob[[j]]
    marker_presence <- as.integer(substr(intersection, j, j))

    if (marker_presence == 1) {
      temp_prob <- temp_prob * marker_prob
    } else {
      temp_prob <- temp_prob * (1 - marker_prob)
    }
  }

  expected_prob_df$Expected_Probability[i] <- temp_prob
}

# Convert expected probabilities to expected counts
expected_prob_df <- expected_prob_df %>%
  mutate(Expected_Count = Expected_Probability * total_observations)

# Observed counts from the actual data
observed <- membership_matrix %>%
  count(Intersection, name = "Observed")

# Simulate null distribution of deviations
set.seed(1009)
null_deviation_distributions <- replicate(1000, {
  # Shuffle the marker expressions
  shuffled <- apply(membership_matrix[, -ncol(membership_matrix)], 2, sample)

  # Create Intersection column from shuffled data
  shuffled_long <- as.data.frame(shuffled) %>%
    mutate(Intersection = apply(., 1, paste0, collapse = ""))

  # Count intersections in shuffled data
  shuffled_counts <- shuffled_long %>%
    count(Intersection, name = "Shuffled_Observed")

  # Join with expected counts to calculate deviations
  deviation <- expected_prob_df %>%
    left_join(shuffled_counts, by = "Intersection") %>%
    mutate(
      Shuffled_Observed = replace_na(Shuffled_Observed, 0),
      Deviation = Shuffled_Observed - Expected_Count
    ) %>%
    select(Intersection, Deviation)

  deviation
}, simplify = FALSE)

# Aggregate null deviations
null_deviation_summary <- bind_rows(null_deviation_distributions, .id = "Simulation") %>%
  group_by(Intersection) %>%
  summarise(
    MeanDeviation = mean(Deviation, na.rm = TRUE),
    SDDeviation = sd(Deviation, na.rm = TRUE),
    .groups = 'drop'
  )

# Calculate deviation statistics for the actual observed data
deviation_stats <- observed %>%
  left_join(expected_prob_df, by = "Intersection") %>%
  mutate(
    Actual_Deviation = Observed - Expected_Count
  ) %>%
  left_join(null_deviation_summary, by = "Intersection") %>%
  mutate(
    p_value = pnorm(Actual_Deviation / SDDeviation, lower.tail = FALSE),
    p_adj = p.adjust(p_value, method = "bonferroni")
  )

# Cap p-values for visualization and avoid -Inf in log transformation
deviation_stats <- deviation_stats %>%
  mutate(
    p_adj_capped = pmax(p_adj, 1e-300),
    log10_p_adj = -log10(p_adj_capped)
  )

# Merge back to membership matrix for plotting
membership_matrix_long <- membership_matrix %>%
  left_join(deviation_stats, by = "Intersection")

# Plot with ComplexUpset
upset_plot <- upset(
  membership_matrix_long,
  intersect = colnames(binary_exp_matrix),
  min_degree = 0,
  annotations = list(
    'Deviation' = (
      ggplot(mapping = aes(fill = Actual_Deviation)) +
        geom_bar() +
        scale_fill_gradient2(
          low = "orange", mid = "pink", high = "darkblue",
          midpoint = 0, na.value = "gray",
          guide = guide_colorbar(title = "Deviation")
        )
    )
  )
)

# Display the plot
print(upset_plot)

deviation_stats_sorted <- deviation_stats[order(-deviation_stats$Observed), ]

library(ggplot2)

# Reorder the 'Intersection' factor based on 'Observed' in descending order
deviation_stats_sorted$Intersection <- factor(deviation_stats_sorted$Intersection,
                                              levels = deviation_stats_sorted$Intersection[order(-deviation_stats_sorted$Observed)])

# Create the plot with the reordered Intersection
ggplot(deviation_stats_sorted, aes(x = Intersection)) +
  theme_classic() +
  geom_point(aes(y = Observed, color = "Observed"), size = 3) +
  geom_point(aes(y = Expected_Count, color = "Expected"), size = 3, shape = 1) +
  scale_color_manual(values = c("Observed" = "pink", "Expected" = "blue")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(#title = "Observed vs Expected Counts",
       y = "Count", color = "Legend")
