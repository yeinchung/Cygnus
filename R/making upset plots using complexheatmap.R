m <- m_list[["Stage 1"]]

UpSet(m[comb_size(m) >= 300])

m <- m[comb_size(m) >= 300]

UpSet(m, top_annotation = HeatmapAnnotation(
  degree = as.character(comb_degree(m)),
  "Intersection\nsize" = anno_barplot(comb_size(m),
                                      border = FALSE,
                                      gp = gpar(fill = "black"),
                                      height = unit(2, "cm")
  ),
  annotation_name_side = "left",
  annotation_name_rot = 0))


combined_data_pre <- combined_data
combined_data <- combined_data[combined_data$`Number.MarkerPositive.Evs.` > 0, ]
states <- combined_data$stage_category
m_list = tapply(seq_len(nrow(combined_data)), states, function(ind) {
  m = make_comb_mat(combined_data[ind, , drop = FALSE])
  m[comb_degree(m) > 0]
  m[comb_size(m) >= 300]
})


sapply(m_list, comb_size)
m_list = normalize_comb_mat(m_list)
sapply(m_list, comb_size)

max_set_size = max(sapply(m_list, set_size))
max_comb_size = max(sapply(m_list, comb_size))

m_list
max_set_size = max(sapply(m_list, set_size))
rel_comb_size = sapply(m_list, function(m) {
  s = comb_size(m)
  s/sum(s)
})



ht_list = NULL
for(i in seq_along(m_list)) {
  ht_list = ht_list %v%
    UpSet(m_list[[i]], row_title = paste0("rating in", names(m_list)[i]),
          set_order = NULL, comb_order = NULL,
          top_annotation = HeatmapAnnotation(
            "Relative\nfraction" = anno_barplot(
              rel_comb_size[, i],
              ylim = c(0, 1),
              gp = gpar(fill = "black"),
              border = FALSE,
              height = unit(2, "cm"),
            ),
            annotation_name_side = "left",
            annotation_name_rot = 0),
          right_annotation = upset_right_annotation(m_list[[i]],
                                                    ylim = c(0, max_set_size))
    )
}
ht_list


###

ht_list = NULL
for(i in seq_along(m_list)) {
  max_rel_fraction = max(rel_comb_size[, i])
  ht_list = ht_list %v%
    UpSet(m_list[[i]], row_title = paste0("rating in ", names(m_list)[i]),
          set_order = NULL, comb_order = NULL,
          top_annotation = HeatmapAnnotation(
            "Relative\nfraction" = anno_barplot(
              rel_comb_size[, i],
              ylim = c(0, max_rel_fraction),
              gp = gpar(fill = "black"),
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






