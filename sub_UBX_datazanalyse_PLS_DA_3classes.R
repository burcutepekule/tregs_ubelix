library(viridis)
library(caret)
library(MASS)
library(mixOmics)
library(car)
library(ggsignif)

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

# ============ DUAL ELLIPSE ANALYSIS ============

# Function to analyze ellipse for a specific class
analyze_class_ellipse = function(plsda_df, target_class, conf_levels) {
  
  # Get points for target class
  class_points = plsda_df %>% filter(treg_outcome == target_class)
  
  # Fit ellipse to class distribution
  class_center = colMeans(class_points[, c("Comp1", "Comp2")])
  class_cov = cov(class_points[, c("Comp1", "Comp2")])
  
  # Calculate Mahalanobis distance for all points
  mahal_dist = mahalanobis(
    plsda_df[, c("Comp1", "Comp2")],
    center = class_center,
    cov = class_cov
  )
  
  ellipse_results = data.frame()
  
  for(conf_level in conf_levels) {
    chi2_threshold = qchisq(conf_level, df = 2)
    in_ellipse = mahal_dist < chi2_threshold
    
    # Calculate metrics
    tp = sum(plsda_df$treg_outcome == target_class & in_ellipse)
    fp = sum(plsda_df$treg_outcome != target_class & in_ellipse)
    fn = sum(plsda_df$treg_outcome == target_class & !in_ellipse)
    
    total_target = sum(plsda_df$treg_outcome == target_class)
    total_other = sum(plsda_df$treg_outcome != target_class)
    
    recall = tp / total_target  # % of target class captured
    precision = tp / (tp + fp)  # % of region that's target class
    pct_other_in = fp / total_other  # % of other classes in region
    enrichment = recall / pct_other_in  # Fold enrichment
    
    ellipse_results = rbind(ellipse_results, data.frame(
      conf_level = conf_level,
      recall = recall,
      precision = precision,
      pct_other_in = pct_other_in,
      enrichment = enrichment,
      n_target_captured = tp,
      n_other_included = fp
    ))
  }
  
  return(list(
    results = ellipse_results,
    center = class_center,
    cov = class_cov,
    mahal_dist = mahal_dist,
    class_points = class_points
  ))
}

# Analyze both classes
conf_levels    = c(0.50, 0.75, 0.90, 0.95, 0.99)

ellipse_plus1  = analyze_class_ellipse(plsda_df, "1", conf_levels)
ellipse_minus1 = analyze_class_ellipse(plsda_df, "-1", conf_levels)

# ============ VISUALIZATIONS ============

# Add ellipse membership to dataframe
plsda_df$in_ellipse_plus1  = ellipse_plus1$mahal_dist < qchisq(level_plus1, df = 2)
plsda_df$in_ellipse_minus1 = ellipse_minus1$mahal_dist < qchisq(level_minus1, df = 2)

# Create a combined category
plsda_df = plsda_df %>%
  mutate(region = case_when(
    in_ellipse_plus1 & in_ellipse_minus1 ~ "Both",
    in_ellipse_plus1 ~ "Tregs better only",
    in_ellipse_minus1 ~ "Tregs worse only",
    TRUE ~ "Tregs don't matter"
  ))

# Function to generate ellipse coordinates
get_ellipse_coords = function(center, cov_matrix, level, n_points = 100) {
  theta = seq(0, 2 * pi, length.out = n_points)
  circle = cbind(cos(theta), sin(theta))
  
  # Scale by chi-square quantile
  radius = sqrt(qchisq(level, df = 2))
  ellipse = circle %*% chol(cov_matrix) * radius
  ellipse = sweep(ellipse, 2, center, "+")
  
  return(data.frame(x = ellipse[, 1], y = ellipse[, 2]))
}

# Generate ellipse coordinates
ellipse_plus1_coords = get_ellipse_coords(
  ellipse_plus1$center, 
  ellipse_plus1$cov, 
  level_plus1
)

ellipse_minus1_coords = get_ellipse_coords(
  ellipse_minus1$center, 
  ellipse_minus1$cov, 
  level_minus1
)

# ============ ADD PLS-DA LOADINGS AS ARROWS ============

# Create a scaling factor for the arrows
arrow_scale = 5  # Adjust this to make arrows longer/shorter

# Get PLS-DA loadings
plsda_loadings = plsda_model$loadings$X

# Create arrows dataframe
plsda_arrows = data.frame(
  parameter = rownames(plsda_loadings),
  Comp1 = plsda_loadings[, 1] * arrow_scale,
  Comp2 = plsda_loadings[, 2] * arrow_scale
)

# Select top N most important parameters
n_arrows = 12 #23 to select all 
plsda_arrows = plsda_arrows %>%
  mutate(total_loading = sqrt(Comp1^2 + Comp2^2)) %>%
  arrange(desc(total_loading)) %>%
  head(n_arrows)

# ============ PLOT WITH ARROWS ============

p2 = ggplot(plsda_df, aes(x = Comp1, y = Comp2)) +
  geom_point(data = plsda_df %>% filter(treg_outcome == "0"),
             aes(color = region), alpha = 0.8, size = 2) +
  scale_color_manual(values = c(
    "Tregs don't matter" = "gray90",
    "Tregs better only" = "lightblue",
    "Tregs worse only" = "pink",
    "Both" = "orange"
  )) +
  # Manual ellipses using exact calculations
  geom_path(data = ellipse_plus1_coords, aes(x = x, y = y),
            color = "blue", linewidth = 1.5, inherit.aes = FALSE) +
  geom_path(data = ellipse_minus1_coords, aes(x = x, y = y),
            color = "red", linewidth = 1.5, inherit.aes = FALSE) +
  geom_point(data = plsda_df %>% filter(treg_outcome == "1"),
             color = "darkblue", alpha = 0.8, size = 3) +
  geom_point(data = plsda_df %>% filter(treg_outcome == "-1"),
             color = "darkred", alpha = 0.8, size = 3) +
  # Add parameter vectors (arrows)
  geom_segment(data = plsda_arrows,
               aes(x = 0, y = 0, xend = Comp1, yend = Comp2),
               arrow = arrow(length = unit(0.3, "cm")),
               color = "black", linewidth = 0.5,
               inherit.aes = FALSE) +
  geom_text_repel(data = plsda_arrows,
                  aes(x = Comp1, y = Comp2, label = parameter),
                  color = "black", size = 3, fontface = "bold",
                  inherit.aes = FALSE) +
  theme_minimal() +
  labs(title = paste0("PLS DA of parameter sets based on Treg help: conf. : (",level_plus1,", ",level_minus1,")"),
       subtitle = sprintf("Enrichment Tregs better: %.1f-fold | Enrichment Tregs worse: %.1f-fold", 
                          ellipse_plus1$results$enrichment[ellipse_plus1$results$conf_level == level_plus1],
                          ellipse_minus1$results$enrichment[ellipse_minus1$results$conf_level == level_minus1]))

print(p2)

ggsave(
  filename = paste0("./treg_regions_PLS_DA_",inj_type,"_w_arrows.png"),
  plot = p2,
  width = 9,
  height = 6,
  dpi = 300,
  bg='white'
)

dev.off()

# ============ PRINT PLS-DA LOADINGS ============

cat("\n========== PLS-DA LOADINGS (Component 1) ==========\n")
plsda_loadings_comp1 = data.frame(
  parameter = rownames(plsda_loadings),
  loading = plsda_loadings[, 1]
) %>%
  arrange(desc(abs(loading)))

print(head(plsda_loadings_comp1, 10))

cat("\n========== PLS-DA LOADINGS (Component 2) ==========\n")
plsda_loadings_comp2 = data.frame(
  parameter = rownames(plsda_loadings),
  loading = plsda_loadings[, 2]
) %>%
  arrange(desc(abs(loading)))

print(head(plsda_loadings_comp2, 10))

# # Plot 3: Separate panels for each class
# 
# p3a = ggplot(plsda_df, aes(x = Comp1, y = Comp2)) +
#   geom_point(data = plsda_df %>% filter(treg_outcome == "0"),
#              aes(color = in_ellipse_plus1), alpha = 0.3, size = 2) +
#   geom_point(data = plsda_df %>% filter(treg_outcome == "1"),
#              color = "darkblue", alpha = 0.8, size = 3) +
#   stat_ellipse(data = ellipse_plus1$class_points, 
#                aes(x = Comp1, y = Comp2),
#                type = "norm", level = level_plus1, 
#                color = "blue", linewidth = 1.5) +
#   scale_color_manual(values = c("gray80", "lightblue")) +
#   theme_minimal() +
#   labs(title = "Class +1 (Tregs Better)",
#        subtitle = sprintf("Enrichment: %.1f-fold", 
#                           ellipse_plus1$results$enrichment[ellipse_plus1$results$conf_level == level_plus1]))
# 
# p3b = ggplot(plsda_df, aes(x = Comp1, y = Comp2)) +
#   geom_point(data = plsda_df %>% filter(treg_outcome == "0"),
#              aes(color = in_ellipse_minus1), alpha = 0.3, size = 2) +
#   geom_point(data = plsda_df %>% filter(treg_outcome == "-1"),
#              color = "darkred", alpha = 0.8, size = 3) +
#   stat_ellipse(data = ellipse_minus1$class_points, 
#                aes(x = Comp1, y = Comp2),
#                type = "norm", level = level_minus1, 
#                color = "red", linewidth = 1.5) +
#   scale_color_manual(values = c("gray80", "pink")) +
#   theme_minimal() +
#   labs(title = "Class -1 (Tregs Worse)",
#        subtitle = sprintf("Enrichment: %.1f-fold", 
#                           ellipse_minus1$results$enrichment[ellipse_minus1$results$conf_level == level_minus1]))
# 
# print(p3a | p3b)
# 
# ggsave(
#   filename = paste0("./treg_regions_PLS_DA_sidebyside.png"),
#   plot = (p3a | p3b),
#   width = 14,
#   height = 6,
#   dpi = 300,
#   bg='white'
# )
# 
# dev.off()

# ============ COMPARISON TABLE ============

# cat("\n========== COMPARISON AT SELECTED CONFIDENCE LEVELS ==========\n")
comparison = data.frame(
  Class = c("+1", "-1"),
  Conf_Level = c(level_plus1, level_minus1),
  N_Total = c(sum(plsda_df$treg_outcome == "1"), 
              sum(plsda_df$treg_outcome == "-1")),
  N_Captured = c(
    ellipse_plus1$results$n_target_captured[ellipse_plus1$results$conf_level == level_plus1],
    ellipse_minus1$results$n_target_captured[ellipse_minus1$results$conf_level == level_minus1]
  ),
  Recall = c(
    ellipse_plus1$results$recall[ellipse_plus1$results$conf_level == level_plus1],
    ellipse_minus1$results$recall[ellipse_minus1$results$conf_level == level_minus1]
  ),
  Precision = c(
    ellipse_plus1$results$precision[ellipse_plus1$results$conf_level == level_plus1],
    ellipse_minus1$results$precision[ellipse_minus1$results$conf_level == level_minus1]
  ),
  Enrichment = c(
    ellipse_plus1$results$enrichment[ellipse_plus1$results$conf_level == level_plus1],
    ellipse_minus1$results$enrichment[ellipse_minus1$results$conf_level == level_minus1]
  )
)

# print(comparison)

# ============ REGION OVERLAP ANALYSIS ============

cat("\n========== REGION OVERLAP ==========\n")
overlap_table = table(
  Plus1_Region = plsda_df$in_ellipse_plus1,
  Minus1_Region = plsda_df$in_ellipse_minus1
)
# print(overlap_table)

# cat("\nPoints in both regions:", sum(plsda_df$in_ellipse_plus1 & plsda_df$in_ellipse_minus1), "\n")
# cat("Points in +1 only:", sum(plsda_df$in_ellipse_plus1 & !plsda_df$in_ellipse_minus1), "\n")
# cat("Points in -1 only:", sum(plsda_df$in_ellipse_minus1 & !plsda_df$in_ellipse_plus1), "\n")
# cat("Points in neither:", sum(!plsda_df$in_ellipse_plus1 & !plsda_df$in_ellipse_minus1), "\n")

# ============ ENRICHMENT COMPARISON PLOTS ============

# # Combine results for plotting
# ellipse_plus1$results$class = "+1"
# ellipse_minus1$results$class = "-1"
# combined_results = rbind(ellipse_plus1$results, ellipse_minus1$results)
# 
# p4 = ggplot(combined_results, aes(x = recall, y = enrichment, color = class)) +
#   geom_line(linewidth = 1.2) +
#   geom_point(size = 3) +
#   scale_color_manual(values = c("+1" = "blue", "-1" = "red")) +
#   theme_minimal() +
#   labs(title = "Enrichment vs Recall for Both Classes",
#        x = "Recall (% of class captured)",
#        y = "Fold Enrichment",
#        color = "Class")
# 
# print(p4)


### REORDER BASED ON WHICH PARAMETER PROFILE YOU WANT (BASED ON COHENS D OR LDA REGION)
### BECAUSE THE LAST SCRIPT WILL DEFINE WHAT "REGION" IS
source("~/Dropbox/tregs_ubelix/sub_UBX_param_hists_cohens_category_3classes.R")
# source("~/Dropbox/tregs_ubelix/sub_UBX_param_hists_DA_region_3classes.R")


# # ======= Step 5: Create Parameter Profiles
# df_profiles_min <- df_analysis %>%
#   dplyr::filter(region %in% c("Tregs better only", "Tregs worse only")) %>%
#   dplyr::select(region, all_of(param_names)) %>%
#   dplyr::group_by(region) %>%
#   summarise(across(all_of(param_names), min)) %>%
#   pivot_longer(cols = -region, names_to = "parameter", values_to = "min_val")
# 
# df_profiles_max <- df_analysis %>%
#   dplyr::filter(region %in% c("Tregs better only", "Tregs worse only")) %>%
#   dplyr::select(region, all_of(param_names)) %>%
#   dplyr::group_by(region) %>%
#   summarise(across(all_of(param_names), max)) %>%
#   pivot_longer(cols = -region, names_to = "parameter", values_to = "max_val")


# ======== Step 6: Interpretable Rules

# # Create decision rules based on regions
# cat("\n========== INTERPRETABLE RULES ==========\n")
# 
# # For "Tregs better" region
# better_rules <- param_tests %>%
#   filter(adj_p_better_vs_worse < 0.05) %>%
#   head(5) %>%
#   mutate(
#     rule = ifelse(mean_better > mean_worse,
#                   paste0(parameter, " HIGH (>", round(mean_better, 2), ")"),
#                   paste0(parameter, " LOW (<", round(mean_better, 2), ")"))
#   )
# 
# cat("\nTregs are BENEFICIAL when:\n")
# cat(paste(better_rules$rule, collapse = "\nAND "))
# 
# # For "Tregs worse" region
# worse_rules <- param_tests %>%
#   filter(adj_p_better_vs_worse < 0.05) %>%
#   head(5) %>%
#   mutate(
#     rule = ifelse(mean_worse > mean_better,
#                   paste0(parameter, " HIGH (>", round(mean_worse, 2), ")"),
#                   paste0(parameter, " LOW (<", round(mean_worse, 2), ")"))
#   )
# 
# cat("\n\nTregs are HARMFUL when:\n")
# cat(paste(worse_rules$rule, collapse = "\nAND "))
# 
# 
# # Create human-readable summary
# summary_table <- param_tests %>%
#   filter(adj_p_better_vs_worse < 0.05) %>%
#   head(10) %>%
#   mutate(
#     better_direction = ifelse(mean_better > mean_worse, "↑ Higher", "↓ Lower"),
#     better_vs_worse = sprintf("%.3f vs %.3f", mean_better, mean_worse),
#     interpretation = sprintf("%s in 'Better' region (%s)", better_direction, better_vs_worse)
#   ) %>%
#   dplyr::select(parameter, interpretation, cohens_d, adj_p_better_vs_worse)
# 
# print("\n========== TOP DIFFERENTIATING PARAMETERS ==========")
# print(summary_table)
# 
# if(dim(summary_table)[1]>1){
#   #=======================
#   # Get only parameters that are significantly different
#   sig_params_ordered <- param_tests %>%
#     dplyr::filter(adj_p_better_vs_worse < 0.05) %>%
#     head(6) %>% # colors depend on this because of z-scoring makes everything relative within each parameter!
#     arrange(adj_p_better_vs_worse) %>%
#     pull(parameter)
#   
#   cat("Significantly different parameters (ordered by p-value):\n")
#   print(sig_params_ordered)
#   
#   # Get profiles for significant parameters only
#   region_profiles_sig <- df_analysis %>%
#     dplyr::filter(region %in% c("Tregs better only", "Tregs worse only", "Tregs don't matter")) %>%
#     dplyr::group_by(region) %>%
#     summarise(across(all_of(sig_params_ordered), mean))
#   
#   profile_matrix_sig <- region_profiles_sig %>%
#     dplyr::select(-region) %>%
#     as.matrix()
#   
#   rownames(profile_matrix_sig) <- region_profiles_sig$region
#   
#   # Transpose and scale (z-score by parameter)
#   profile_matrix_sig_scaled <- scale(t(profile_matrix_sig))
#   
#   # CRITICAL: Reorder rows to match sig_params_ordered (no clustering!)
#   profile_matrix_sig_scaled <- profile_matrix_sig_scaled[sig_params_ordered, ]
#   
#   # Save the heatmap
#   png(filename = paste0("./treg_params_heatmap_PLS_DA_",inj_type,".png"),
#       width = 6,
#       height = 4,
#       units = "in",
#       res = 300,
#       bg = "white")
#   
#   pheatmap(profile_matrix_sig_scaled,
#            cluster_rows = FALSE,
#            cluster_cols = FALSE,
#            color = colorRampPalette(c("blue", "white", "red"))(100),
#            main = "Significantly Different Parameters\n(ordered by p-value, z-scored)",
#            fontsize = 10,
#            cellwidth = 80,
#            cellheight = 15)
#   
#   dev.off()
# }



