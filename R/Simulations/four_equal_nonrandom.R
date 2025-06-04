#install.packages("ggplot2")
#install.packages("ggthemes")
#install.packages("tidyverse")
#install.packages("ComplexUpset")

library(ggplot2)
library(ggrepel)
library(ggthemes)
library(scales)
library(tidyverse)
library(ComplexUpset)

# code for four markers, equal proportions, A and B more likely to co-localize
set.seed(42)

columns <- c('A', 'B', 'C', 'D')
proportions <- c(0.2, 0.2, 0.2, 0.2)

rows <- 4000

# ======= make binary matrix =======

binary_exp_matrix <- data.frame(matrix(ncol = length(columns), nrow = rows))

make_columns <- function(proportion, rows) {
  num_ones <- floor(rows * proportion)
  column_data <- c( rep(1, num_ones), rep(0, rows - num_ones))
  return(sample(column_data))
}

for (i in 1:length(columns)) {
  binary_exp_matrix[, i] <- make_columns(proportions[i], rows)
}

# adjust probabilities
# shift_factor <-0.8 # this can be adjusted and later function can be set to runif(1) < shift_factor
# to change the "degree of correlation"...

colnames(binary_exp_matrix) <- columns
A_indices <- which(binary_exp_matrix$A == 1)

# this version makes all rows where A = 1 also have B = 1
for (idx in A_indices) {
  if (TRUE) { #runif(1) < shift_factor
    binary_exp_matrix[idx, "B"] <- 1
    binary_exp_matrix[idx, "C"] <- 1
    binary_exp_matrix[idx, "D"] <- 1
  }
}

colSums(binary_exp_matrix) # can note that the proportions are messed up...

# this fixes proportions randomly
for (col in c("B", "C", "D")) {
  ones_to_remove <- sum(binary_exp_matrix[, col]) - floor(rows * proportions[which(columns == col)])
  if (ones_to_remove > 0) {
    ones_positions <- which(binary_exp_matrix[, col] == 1)
    remove_indices <- sample(ones_positions, ones_to_remove)
    binary_exp_matrix[remove_indices, col] <- 0
  }
}

colSums(binary_exp_matrix) # proportions should be fixed!

# add a column to add intersections
membership_matrix <- as.data.frame(binary_exp_matrix) %>%
  mutate(Intersection = apply(., 1, paste0, collapse = ""))  # Create a unique identifier for each combination of markers


total_obsv <- nrow(membership_matrix) # this is same as rows previously set !


# calculate the probability of each marker being expressed (1)
prob <- membership_matrix %>%
  summarise(across(everything(), ~ mean(. == 1)  ))

within(prob, rm(Intersection))


# ======= make new df with expected probabilities  =======
present_intersections <- unique(membership_matrix$Intersection)
expected_prob_df <- data.frame(Intersection = present_intersections, Expected_Probability = NA)


prob <- within(prob, rm(Intersection))

# calculate probability for each intersection
for (i in seq_along(present_intersections)) {
  intersection <- present_intersections[i]
  temp_prob <- 1

  for (j in seq_along(prob)) {
    marker_prob <- prob[[j]]  # probability of the j-th marker being expressed
    marker_presence <- as.integer(substr(intersection, j, j))  # presence (1) or absence (0) of the marker in the intersection

    if (marker_presence == 1) {
      temp_prob <- temp_prob * marker_prob   # multiply by the probability of presence
    } else {
      temp_prob  <- temp_prob * (1 - marker_prob)  # m ltiply by the probability of absence
    }
  }
  expected_prob_df$Expected_Probability[i] <- temp_prob
}


# get expected counts
expected_prob_df <- expected_prob_df %>%
  mutate(Expected_Count = Expected_Probability * total_obsv)

# get observed counts
observed <- membership_matrix %>%
  count(Intersection, name = "Observed")

# ======= simulation =======

# simulate null distribution of deviations by shuffling marker expressions
iteration_no = 1000

set.seed(04)
null_deviation_distributions <- replicate(iteration_no, {
  # simulate random co-expression
  shuffled <- apply(membership_matrix[, -ncol(membership_matrix)], 2, sample)

  shuffled_long <- as.data.frame(shuffled) %>%
    mutate(Intersection = apply(., 1, paste0, collapse = ""))

  shuffled_counts <- shuffled_long %>%
    count(Intersection, name = "Shuffled_Observed")

  deviation <- expected_prob_df %>%
    left_join(shuffled_counts, by = "Intersection") %>%
    mutate(
      Shuffled_Observed = replace_na(Shuffled_Observed, 0),
      Deviation = Shuffled_Observed - Expected_Count
    ) %>%
    select(Intersection, Deviation)

  deviation
}, simplify = FALSE)


# aggregate null deviations across simulations
null_deviation_summary <- bind_rows(null_deviation_distributions, .id = "Simulation") %>%
  group_by(Intersection) %>%
  summarise(
    MeanDeviation = mean(Deviation, na.rm = TRUE),  # mean deviation
    SDDeviation = sd(Deviation, na.rm = TRUE),      # sd deviation
    .groups = 'drop'
  )

# calculate deviation statistics
deviation_stats <- observed %>%
  left_join(expected_prob_df, by = "Intersection") %>%
  mutate(
    Actual_Deviation = Observed - Expected_Count    # compute actual deviation
  ) %>%
  left_join(null_deviation_summary, by = "Intersection") %>%
  mutate(
    p_value = pnorm(Actual_Deviation / SDDeviation, lower.tail = FALSE),   # compute p-value
    p_adj = p.adjust(p_value, method = "bonferroni")   # important for multiple testing and false pos!
  )


# cap p-val to avoid errors later
deviation_stats <- deviation_stats %>%
  mutate(
    p_adj_capped = pmax(p_adj, 1e-300),
    log10_p_adj = -log10(p_adj_capped)     # this for volcano plot
  )

membership_matrix_long <- membership_matrix %>%
  left_join(deviation_stats, by = "Intersection")


# ======= make plots ! =======
membership_matrix_mixed <- membership_matrix_long %>%
  mutate(
    significant = ifelse(is.na(p_adj), "",  # add significance labels
                         ifelse(p_adj < 0.005, "*", ""))
  )


upset_plot <-upset(
  membership_matrix_mixed,
  intersect = colnames(binary_exp_matrix),
  min_degree=1,
  name="biomarkers",
  stripes = c('cornsilk1', 'deepskyblue1', 'grey90'),
  set_sizes = (
    upset_set_size()
    + geom_text(
      aes(label = ..count..),
      stat = "count",
      hjust = -0.5,
      size = 3,
      color = "white"
    )
    + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  ),

  base_annotations=list(
    'Intersection size'=intersection_size(
      counts= TRUE,
      bar_number_threshold = 1,
      mapping =aes(fill = Actual_Deviation),
      text_colors=c(
        on_background='brown', on_bar='yellow'
      ),
    )
    + annotate(
      geom='text', x=Inf, y=Inf,
      label = paste('* : p-value < 0.005' ),
      vjust=1.5, hjust=1
    )
    + geom_text(
      aes(
        label = significant,
        y = Observed +50
      ),
      size = 5,
      vjust = -0.5,
      color = "black"
    )
    + scale_fill_gradient2(
      low = "orange", mid = "pink", high = "darkblue",
      midpoint = 0, na.value = "gray",
      guide = guide_colorbar(title = "Deviation")
    )
    + ylab('Intersection size')
    + scale_y_continuous(expand = expansion(mult = c(0, 0.2)))
  ),

  width_ratio = 0.1
)

print(upset_plot)

# volcano
significance_threshold <- 0.05
log_threshold <- -log10(significance_threshold)

ggplot(deviation_stats, aes(x = Actual_Deviation, y = log10_p_adj, label = Intersection)) +
  geom_point(aes(color = log10_p_adj, size = abs(Actual_Deviation)), alpha = 0.7) +
  geom_text_repel(data = subset(deviation_stats, log10_p_adj > log_threshold),
                  aes(label = Intersection),
                  size = 5, fontface = "bold", color = "darkblue",
                  box.padding = 0.6, max.overlaps = 15) +
  geom_hline(yintercept = log_threshold, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_gradientn(colors = c("blue", "purple", "red"),
                        limits = c(0, max(deviation_stats$log10_p_adj)), oob = squish) +
  labs(
    x = "Actual Deviation",
    y = expression(-log[10](adjusted~p~value)),
    color = "-log10(p_adj)",
    size = "Deviation Magnitude"
  ) +
  theme_classic(base_size = 15) +
  theme(legend.position = "right")
