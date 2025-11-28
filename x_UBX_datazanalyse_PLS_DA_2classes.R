plsda_model = plsda(
  X = df_lda %>% dplyr::select(all_of(param_names)),
  Y = df_lda$tregs_better_cohens,
  ncomp = 2
)

# Get PLS-DA projections
plsda_proj = plsda_model$variates$X
plsda_df = data.frame(
  Comp1 = plsda_proj[, 1],
  Comp2 = plsda_proj[, 2],
  treg_outcome = as.factor(df_lda$tregs_better_cohens)
)

# ============ ELLIPSE ANALYSIS ============

# Get class 1 points
class1_points = plsda_df %>% filter(treg_outcome == "1")

# Fit ellipse to class 1 distribution
class1_center = colMeans(class1_points[, c("Comp1", "Comp2")])
class1_cov = cov(class1_points[, c("Comp1", "Comp2")])

# Calculate Mahalanobis distance for all points
mahal_dist = mahalanobis(
  plsda_df[, c("Comp1", "Comp2")],
  center = class1_center,
  cov = class1_cov
)

# Test different confidence levels and calculate enrichment
cat("\n========== ELLIPSE ENRICHMENT ANALYSIS ==========\n")

ellipse_results = data.frame()

for(conf_level in c(0.50, 0.75, 0.90, 0.95, 0.99)) {
  chi2_threshold = qchisq(conf_level, df = 2)
  in_ellipse = mahal_dist < chi2_threshold
  
  # Calculate metrics
  tp = sum(plsda_df$treg_outcome == "1" & in_ellipse)
  fp = sum(plsda_df$treg_outcome == "0" & in_ellipse)
  fn = sum(plsda_df$treg_outcome == "1" & !in_ellipse)
  tn = sum(plsda_df$treg_outcome == "0" & !in_ellipse)
  
  total_1s = sum(plsda_df$treg_outcome == "1")
  total_0s = sum(plsda_df$treg_outcome == "0")
  
  recall = tp / total_1s  # % of 1s captured
  precision = tp / (tp + fp)  # % of region that's 1s
  pct_0s_in = fp / total_0s  # % of 0s in region
  enrichment = recall / pct_0s_in  # Fold enrichment
  
  ellipse_results = rbind(ellipse_results, data.frame(
    conf_level = conf_level,
    recall = recall,
    precision = precision,
    pct_0s_in = pct_0s_in,
    enrichment = enrichment,
    n_1s_captured = tp,
    n_0s_included = fp
  ))
  
  cat(sprintf("\n%d%% Confidence Ellipse:\n", conf_level * 100))
  cat(sprintf("  Captures: %d/%d class 1 points (%.1f%% recall)\n", 
              tp, total_1s, recall * 100))
  cat(sprintf("  Includes: %d/%d class 0 points (%.1f%%)\n", 
              fp, total_0s, pct_0s_in * 100))
  cat(sprintf("  Precision: %.1f%%\n", precision * 100))
  cat(sprintf("  Enrichment: %.1f-fold\n", enrichment))
}

print("\nSummary Table:")
print(ellipse_results)

# ============ VISUALIZATIONS ============

# Choose 100*level_in% confidence level for visualization
level_in = 0.75
plsda_df$in_ellipse = mahal_dist < qchisq(level_in, df = 2)

# Plot 1: Basic ellipse visualization
p1 = ggplot(plsda_df, aes(x = Comp1, y = Comp2)) +
  stat_ellipse(data = class1_points, 
               aes(x = Comp1, y = Comp2),
               type = "norm", level = level_in, 
               color = "red", linewidth = 1.2) +
  geom_point(data = plsda_df %>% filter(treg_outcome == "0"),
             aes(color = in_ellipse), alpha = 0.3, size = 2) +
  geom_point(data = plsda_df %>% filter(treg_outcome == "1"),
             color = "blue", alpha = 0.8, size = 3) +
  scale_color_manual(values = c("gray80", "orange"),
                     name = "Class 0",
                     labels = c("Outside", "Inside")) +
  theme_minimal() +
  labs(title = paste0(100*level_in, "% Confidence Ellipse Around Class 1"),
       subtitle = sprintf("Enrichment: %.1f-fold", 
                          ellipse_results$enrichment[ellipse_results$conf_level == level_in]))

print(p1)

# # Plot 2: Multiple ellipse levels
# p2 = ggplot(plsda_df, aes(x = Comp1, y = Comp2)) +
#   geom_point(data = plsda_df %>% filter(treg_outcome == "0"),
#              color = "gray80", alpha = 0.2, size = 2) +
#   stat_ellipse(data = class1_points, 
#                aes(x = Comp1, y = Comp2),
#                type = "norm", level = 0.50, 
#                color = "#CB0404", linewidth = 1, linetype = "solid") +
#   stat_ellipse(data = class1_points, 
#                aes(x = Comp1, y = Comp2),
#                type = "norm", level = 0.75, 
#                color = "#F4631E", linewidth = 1, linetype = "solid") +
#   stat_ellipse(data = class1_points, 
#                aes(x = Comp1, y = Comp2),
#                type = "norm", level = 0.95, 
#                color = "#FF9F00", linewidth = 1, linetype = "solid") +
#   geom_point(data = plsda_df %>% filter(treg_outcome == "1"),
#              color = "#003092", alpha = 0.8, size = 3) +
#   theme_minimal() +
#   labs(title = "Multiple Ellipse Confidence Levels",
#        subtitle = "Red = 50%, Orange = 75%, Blue = 95%")
# 
# print(p2)
# 
# # Plot 3: Enrichment metrics
# 
# p3a = ggplot(ellipse_results, aes(x = recall, y = precision)) +
#   geom_line(linewidth = 1.2) +
#   geom_point(size = 3) +
#   geom_text(aes(label = sprintf("%d%%", conf_level * 100)),
#             nudge_y = 0.01, size = 3) +
#   theme_minimal() +
#   labs(title = "Precision-Recall Trade-off",
#        x = "Recall (% of 1s captured)",
#        y = "Precision (% of region that's 1s)")
# 
# p3b = ggplot(ellipse_results, aes(x = recall, y = enrichment)) +
#   geom_line(linewidth = 1.2, color = "darkblue") +
#   geom_point(size = 3, color = "darkblue") +
#   geom_text(aes(label = sprintf("%d%%", conf_level * 100)),
#             nudge_y = 1, size = 3) +
#   theme_minimal() +
#   labs(title = "Enrichment vs Recall",
#        x = "Recall (% of 1s captured)",
#        y = "Fold Enrichment")
# 
# print(p3a / p3b)
# 
# # ============ FINAL SUMMARY ============
# 
# best_enrichment = ellipse_results[which.max(ellipse_results$enrichment), ]
# cat("\n========== BEST ENRICHMENT ==========\n")
# cat(sprintf("Best enrichment at %d%% confidence level\n", best_enrichment$conf_level * 100))
# cat(sprintf("  %.1f-fold enrichment\n", best_enrichment$enrichment))
# cat(sprintf("  Captures %.1f%% of class 1\n", best_enrichment$recall * 100))
# cat(sprintf("  Precision: %.1f%%\n", best_enrichment$precision * 100))
# 
