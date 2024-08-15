# Function to calculate deviation for each exclusive intersection
calc_deviation <- function(comb_matrix, total_size) {
  deviations <- numeric(ncol(comb_matrix))
  for (i in seq_len(ncol(comb_matrix))) {
    I <- sum(comb_matrix[, i])  # Number of elements in the intersection I
    S_plus <- which(comb_matrix[, i] == 1)  # Sets contained in the exclusive intersection
    S_minus <- which(comb_matrix[, i] == 0)  # Sets not contained in the intersection

    # Calculate products for sets in S*+ and S*-
    if (length(S_plus) > 0) {
      prod_plus <- prod(sapply(S_plus, function(j) sum(comb_matrix[, j]) / total_size))
    } else {
      prod_plus <- 1
    }

    if (length(S_minus) > 0) {
      prod_minus <- prod(sapply(S_minus, function(j) 1 - sum(comb_matrix[, j]) / total_size))
    } else {
      prod_minus <- 1
    }

    # Calculate expected value and deviation
    expected_I <- total_size * prod_plus * prod_minus
    deviations[i] <- abs(I / total_size - expected_I)
  }
  return(deviations)
}





# Calculate deviations for each m in m_list
deviation_list <- lapply(m_list, function(m) {
  comb_matrix <- as.matrix(m)
  total_size <- nrow(comb_matrix)
  calc_deviation(comb_matrix, total_size)
})

# Define a color function for deviation
library(circlize)
col_fun <- colorRamp2(c(min(unlist(deviation_list)), max(unlist(deviation_list))), c("blue", "red"))

ht_list <- NULL
for (i in seq_along(m_list)) {
  max_rel_fraction <- max(rel_comb_size[, i])
  deviation_colors <- col_fun(deviation_list[[i]])
  ht_list <- ht_list %v%
    UpSet(m_list[[i]], row_title = paste0("rating in ", names(m_list)[i]),
          set_order = NULL, comb_order = NULL,
          top_annotation = HeatmapAnnotation(
            "Relative\nfraction" = anno_barplot(
              rel_comb_size[, i],
              ylim = c(0, max_rel_fraction),
              gp = gpar(fill = deviation_colors),
              border = FALSE,
              height = unit(2, "cm")
            ),
            annotation_name_side = "left",
            annotation_name_rot = 0),
          right_annotation = upset_right_annotation(m_list[[i]],
                                                    ylim = c(0, max_set_size))
    )
}
ht_list

