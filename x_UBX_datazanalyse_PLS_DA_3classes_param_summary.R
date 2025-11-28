source("~/Dropbox/tregs_ubelix/UBX_datazanalyse_PLS_DA_3classes.R")

# Create human-readable summary
summary_table <- param_tests %>%
  filter(adj_p_better_vs_worse < 0.05) %>%
  head(10) %>%
  mutate(
    better_direction = ifelse(mean_better > mean_worse, "↑ Higher", "↓ Lower"),
    better_vs_worse = sprintf("%.3f vs %.3f", mean_better, mean_worse),
    interpretation = sprintf("%s in 'Better' region (%s)", better_direction, better_vs_worse)
  ) %>%
  dplyr::select(parameter, interpretation, cohens_d, adj_p_better_vs_worse)

print("\n========== TOP DIFFERENTIATING PARAMETERS ==========")
print(summary_table)

#=======================
# Get only parameters that are significantly different
sig_params_ordered <- param_tests %>%
  dplyr::filter(adj_p_better_vs_worse < 0.05) %>%
  head(10) %>% # colors depend on this because of z-scoring makes everything relative within each parameter!
  arrange(adj_p_better_vs_worse) %>%
  pull(parameter)

cat("Significantly different parameters (ordered by p-value):\n")
print(sig_params_ordered)

# Get profiles for significant parameters only
region_profiles_sig <- df_analysis %>%
  dplyr::filter(region %in% c("Tregs better only", "Tregs worse only", "Tregs don't matter")) %>%
  dplyr::group_by(region) %>%
  summarise(across(all_of(sig_params_ordered), mean))

profile_matrix_sig <- region_profiles_sig %>%
  dplyr::select(-region) %>%
  as.matrix()

rownames(profile_matrix_sig) <- region_profiles_sig$region

# Transpose and scale (z-score by parameter)
profile_matrix_sig_scaled <- scale(t(profile_matrix_sig))

# CRITICAL: Reorder rows to match sig_params_ordered (no clustering!)
profile_matrix_sig_scaled <- profile_matrix_sig_scaled[sig_params_ordered, ]

library(pheatmap)
pheatmap(profile_matrix_sig_scaled,
         cluster_rows = FALSE,  # Turn OFF clustering
         cluster_cols = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Significantly Different Parameters\n(ordered by p-value, z-scored)",
         fontsize = 11,
         cellwidth = 100,
         cellheight = 20)