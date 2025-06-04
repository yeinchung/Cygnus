datadir <- "/Users/yeinchung/Documents/hmb265/yein_mouse/csvs/"
# Load necessary libraries
library(ggplot2)      # For data visualization
library(dplyr)        # For data manipulation
library(tidyverse)    # For data wrangling and visualization
library(ComplexUpset) # For creating complex upset plots
library(Cygnus)       # For single-cell data analysis

setwd(datadir)
mouse1 <- read.csv("Mouse1_Want_20230823.csv")
mouse2 <- read.csv("Mouse2_Want_20230823.csv")
mouse4 <- read.csv("Mouse4_Want_20230823.csv")

mouse1[["Mouse"]] <- 1
mouse2[["Mouse"]] <- 2
mouse4[["Mouse"]] <- 4

merged <- rbind(mouse1, mouse2)
merged <- rbind(merged, mouse4)
dim(merged)
dim(mouse1)

merged <- merged[sample(nrow(merged), 5000), ]

write.csv(merged, "merged.csv")

# Define the markers of interest
raw_markers <-  colnames(merged)[8:17]
meta <- colnames(merged)[18:20]

# Create a Cygnus object using the data file and specified markers
cyg <- CreateCygnus("merged.csv",
                    markers_col = raw_markers,
                    meta_col = meta)

threshold <- c(111, 43, 113, 112, 43, 144, 102, 28, 797, 70)
for(i in 1:length(threshold)){
  threshold[i] = threshold[i] * 1.5
}
threshold

# Extract a binary expression matrix based on a threshold of 1000
#cyg <- scaleExpressionMatrix(cyg)
#cyg <- Cygnus::markRelevantMarkersUI(cyg)
cyg <- Cygnus::createBinaryMatrix(cyg, thresholds = threshold)
#binary_exp_matrix <- cyg@matrices[["binary_exp_matrix"]][, c("EpCAM", "EGFR", "SDC1", "MET", "ADAM10")]
binary_exp_matrix <- cyg@matrices[["binary_exp_matrix"]]

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
iteration_no = 1000 # Number of simulations (1000 for final analysis)

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


# Generate upset plot
upset_plot <-upset(
  membership_matrix_mixed,
  intersect = colnames(binary_exp_matrix),   # Specify markers to include in the plot
  min_degree=2,              # Include all intersections
  name="biomarkers",         # Name for the plot
  stripes = c('cornsilk1', 'deepskyblue1', 'grey90'),  # Color stripes for the plot
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
          y = Observed +50
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
    + scale_y_continuous(expand = expansion(mult = c(0, 0.2)))   # Adjust y-axis scale
  ),

  width_ratio = 0.1,  # Adjust width ratio of the plot
  min_size = 15

)



# Display the upset plot
print(upset_plot)









