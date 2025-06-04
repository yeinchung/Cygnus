# Load required libraries
library(ComplexUpset)
library(ggplot2)
library(dplyr)
library(tidyr)

# Extract the binary expression matrix (excluding the first column)
binary_exp_matrix <- cyg@matrices[["binary_exp_matrix"]][, c("EpCAM", "EGFR", "SDC1", "MET", "ADAM10")]

# Convert the binary matrix into a format suitable for UpSet
binary_long <- binary_exp_matrix %>%
  mutate(Element = rownames(.)) %>%
  pivot_longer(-Element, names_to = "Set", values_to = "Membership") %>%
  pivot_wider(names_from = "Set", values_from = "Membership")

# Compute observed cardinality for each intersection
observed <- binary_exp_matrix %>%
  rowwise() %>%
  mutate(Intersection = paste0(c_across(), collapse = "")) %>%
  ungroup() %>%
  count(Intersection, name = "Observed")

# Simulate null distributions for expected cardinality
set.seed(42)
null_distributions <- replicate(1000, {
  shuffled <- apply(binary_exp_matrix, 2, sample)
  shuffled_intersections <- apply(shuffled, 1, paste0, collapse = "")
  as.data.frame(table(shuffled_intersections)) %>%
    rename(Intersection = shuffled_intersections, Expected = Freq)
}, simplify = FALSE)

# Summarize null distributions
null_distribution_summary <- do.call(bind_rows, null_distributions) %>%
  group_by(Intersection) %>%
  summarise(
    MeanExpected = mean(Expected, na.rm = TRUE),
    SDExpected = sd(Expected, na.rm = TRUE)
  )

# Merge observed and null distribution data
# Prevent -log10(0) by adding a small constant to p-values
deviation_stats <- deviation_stats %>%
  mutate(
    p_adj_capped = pmax(p_adj, 1e-300),  # Cap p_adj to avoid 0
    log10_p_adj = -log10(p_adj_capped)  # Compute -log10(p_adj)
  )

# Set a maximum value for log10(p_adj) to avoid Inf in visualization
max_log10 <- 300
deviation_stats <- deviation_stats %>%
  mutate(log10_p_adj = pmin(log10_p_adj, max_log10))  # Cap at 300 for visualization

# Merge back to membership_matrix
membership_matrix$Intersection <- membership_matrix_long$Intersection
membership_matrix <- membership_matrix %>%
  left_join(deviation_stats, by = "Intersection")


# Create UpSet plot with adjusted -log10(p_adj) for the subset
upset_plot <- upset(
  binary_exp_matrix,
  intersect = colnames(binary_exp_matrix)[!colnames(binary_exp_matrix) %in% "Intersection"],
  annotations = list(
    'Deviation' = (
      ggplot(deviation_stats_subset, aes(fill = deviation_stats_subset$log10_p_adj)) +
        geom_bar(stat = 'count') +
        scale_fill_gradient2(
          low = "blue", mid = "yellow", high = "red",
          midpoint = -log10(0.05), na.value = "gray",
          limits = c(0, max_log10),  # Adjust scale limits as needed
          guide = guide_colorbar(title = "-log10(p-value)")
        )
    )
  )
)

# Display the plot
print(upset_plot)
