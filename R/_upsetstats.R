# libraries
library(ggplot2)
library(dplyr)
library(tidyverse)

# this extracts binary exp matrix from the example Cygnus object
cyg <- Cygnus::createBinaryMatrix(cyg, thresholds = 1000)
binary_exp_matrix <- cyg@matrices[["binary_exp_matrix"]][, c("EpCAM", "EGFR", "SDC1", "MET", "ADAM10")]

# Compute observed cardinality
membership_matrix <- as.data.frame(binary_exp_matrix)
    # > head(membership_matrix)
    # EpCAM EGFR SDC1 MET ADAM10
    # 1     0    0    0   0      0
    # 2     0    0    0   0      0
    # 3     0    0    0   0      0
    # 4     0    0    0   0      0
    # 5     0    0    0   0      0
    # 6     0    0    0   0      0
membership_matrix_long <- membership_matrix %>%
  mutate(Intersection = apply(membership_matrix, 1, paste0, collapse = ""))
    # > head(membership_matrix_long)
    # EpCAM EGFR SDC1 MET ADAM10 Intersection
    # 1     0    0    0   0      0        00000
    # 2     0    0    0   0      0        00000
    # 3     0    0    0   0      0        00000
    # 4     0    0    0   0      0        00000
    # 5     0    0    0   0      0        00000
    # 6     0    0    0   0      0        00000

observed <- membership_matrix_long %>%
  group_by(Intersection) %>%
  summarise(Observed = n())
    # > observed
    # # A tibble: 11 × 2
    # Intersection Observed
    # <chr>           <int>
    #   1 00000            1928
    # 2 00100               3
    # 3 00111               1
    # 4 01111             151
    # 5 10000             358
    # 6 10001               1
    # 7 10100               2
    # 8 10111               5
    # 9 11011               1
    # 10 11100               1
    # 11 11111             813

# simulate null distribution of expected cardinalities
set.seed(1009)
null_distributions <- replicate(1000, {
  shuffled <- apply(membership_matrix, 2, sample)
  shuffled_long <- as.data.frame(shuffled) %>%
    mutate(Intersection = apply(shuffled, 1, paste0, collapse = ""))
  shuffled_long %>%
    group_by(Intersection) %>%
    summarise(Expected = n())
}, simplify = FALSE)
    # > head(null_distributions)
    # [[1]]
    # # A tibble: 32 × 2
    # Intersection Expected
    # <chr>           <int>
    #   1 00000             542
    # 2 00001             206
    # 3 00010             196
    # 4 00011              90
    # 5 00100             186
    # 6 00101             106
    # 7 00110              89
    # 8 00111              45
    # 9 01000             205
    # 10 01001              88
    # # ℹ 22 more rows
    # # ℹ Use `print(n = ...)` to see more rows

# this calculate deviation and p-values
null_distribution_summary <- do.call(bind_rows, null_distributions) %>%
  group_by(Intersection) %>%
  summarise(
    MeanExpected = mean(Expected, na.rm = TRUE),
    SDExpected = sd(Expected, na.rm = TRUE)
  )
    # # A tibble: 32 × 3
    # Intersection MeanExpected SDExpected
    # <chr>               <dbl>      <dbl>
    #   1 00000               507.       16.0
    # 2 00001               216.       11.8
    # 3 00010               215.       12.0
    # 4 00011                91.1       8.42
    # 5 00100               217.       11.8
    # 6 00101                91.5       8.28
    # 7 00110                91.9       8.45
    # 8 00111                38.9       6.00
    # 9 01000               214.       11.7
    # 10 01001                90.0       8.69
    # # ℹ 22 more rows
    # # ℹ Use `print(n = ...)` to see more rows

deviation_stats <- observed %>%
  left_join(null_distribution_summary, by = "Intersection") %>%
  mutate(
    Deviation = Observed - MeanExpected,
    p_value = pnorm(Deviation / SDExpected, lower.tail = FALSE)
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "bonferroni"))
    # > head(deviation_stats)
    # # A tibble: 6 × 7
    # Intersection Observed MeanExpected SDExpected Deviation   p_value     p_adj
    # <chr>           <int>        <dbl>      <dbl>     <dbl>     <dbl>     <dbl>
    #   1 00000            1928        507.       16.0     1421.  0         0
    # 2 00100               3        217.       11.8     -214.  1   e+  0 1   e+  0
    # 3 00111               1         38.9       6.00     -37.9 1.00e+  0 1   e+  0
    # 4 01111             151         16.5       3.85     134.  5.14e-268 5.66e-267
    # 5 10000             358        288.       12.7       70.4 1.42e-  8 1.56e-  7
    # 6 10001               1        122.        9.24    -121.  1   e+  0 1   e+  0

# NEED TO EDITL this caps p-value and log10 accordingly
deviation_stats <- deviation_stats %>%
  mutate(
    p_adj_capped = pmax(p_adj, 1e-300),
    log10_p_adj = -log10(p_adj_capped)
  )

# merge back to membership matrix
membership_matrix_long <- membership_matrix_long %>%
  left_join(deviation_stats, by = "Intersection")

# plot!
upset_plot <- upset(
  membership_matrix_long,
  intersect = colnames(binary_exp_matrix),min_degree=1,
  annotations = list(
    'Deviation' = (
      ggplot(mapping = aes(fill = Deviation)) + # theme_classic() + #log10_p_adj
        geom_bar() + #stat = 'count'
        scale_fill_gradient2(
           low="orange", mid="pink",high = "blue",
          midpoint = 0, na.value = "gray",
          guide = guide_colorbar(title = "Deviation")
        )
    )
  )
)

# display
print(upset_plot)

