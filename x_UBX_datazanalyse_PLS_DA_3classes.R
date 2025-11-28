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
  
  cat(sprintf("\n========== ELLIPSE ANALYSIS FOR CLASS %s ==========\n", target_class))
  
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
    
    cat(sprintf("\n%d%% Confidence Ellipse:\n", conf_level * 100))
    cat(sprintf("  Captures: %d/%d class %s points (%.1f%% recall)\n", 
                tp, total_target, target_class, recall * 100))
    cat(sprintf("  Includes: %d/%d other class points (%.1f%%)\n", 
                fp, total_other, pct_other_in * 100))
    cat(sprintf("  Precision: %.1f%%\n", precision * 100))
    cat(sprintf("  Enrichment: %.1f-fold\n", enrichment))
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
conf_levels = c(0.50, 0.75, 0.90, 0.95, 0.99)

ellipse_plus1 = analyze_class_ellipse(plsda_df, "1", conf_levels)
ellipse_minus1 = analyze_class_ellipse(plsda_df, "-1", conf_levels)

# ============ VISUALIZATIONS ============

# Choose confidence levels for visualization
level_plus1  = 0.75
level_minus1 = 0.75

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

if(dim(ellipse_minus1$class_points)[1]>0){
  ellipse_minus1_coords = get_ellipse_coords(
    ellipse_minus1$center, 
    ellipse_minus1$cov, 
    level_minus1
  )
  
  # Then in your plot, replace stat_ellipse with geom_path:
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
    theme_minimal() +
    labs(title = "PLS DA of parameter sets based on Treg help.",
         subtitle = sprintf("Enrichment Tregs better: %.1f-fold | Enrichment Tregs worse: %.1f-fold", 
                            ellipse_plus1$results$enrichment[ellipse_plus1$results$conf_level == level_plus1],
                            ellipse_minus1$results$enrichment[ellipse_minus1$results$conf_level == level_minus1]))
}else{
  # Then in your plot, replace stat_ellipse with geom_path:
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
    geom_point(data = plsda_df %>% filter(treg_outcome == "1"),
               color = "darkblue", alpha = 0.8, size = 3) +
    theme_minimal() +
    labs(title = "PLS DA of parameter sets based on Treg help.",
         subtitle = sprintf("Enrichment Tregs better: %.1f-fold", 
                            ellipse_plus1$results$enrichment[ellipse_plus1$results$conf_level == level_plus1]))
}


print(p2)

ggsave(
  filename = paste0("./treg_regions_PLS_DA.png"),
  plot = p2,
  width = 9,
  height = 6,
  dpi = 300,
  bg='white'
)

dev.off()

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

cat("\n========== COMPARISON AT SELECTED CONFIDENCE LEVELS ==========\n")
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

print(comparison)

# ============ REGION OVERLAP ANALYSIS ============

cat("\n========== REGION OVERLAP ==========\n")
overlap_table = table(
  Plus1_Region = plsda_df$in_ellipse_plus1,
  Minus1_Region = plsda_df$in_ellipse_minus1
)
print(overlap_table)

cat("\nPoints in both regions:", sum(plsda_df$in_ellipse_plus1 & plsda_df$in_ellipse_minus1), "\n")
cat("Points in +1 only:", sum(plsda_df$in_ellipse_plus1 & !plsda_df$in_ellipse_minus1), "\n")
cat("Points in -1 only:", sum(plsda_df$in_ellipse_minus1 & !plsda_df$in_ellipse_plus1), "\n")
cat("Points in neither:", sum(!plsda_df$in_ellipse_plus1 & !plsda_df$in_ellipse_minus1), "\n")

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

# ========= Step 1: Extract Parameter Values for Each Region
# Add region membership and original parameters to analysis dataframe
df_analysis <- df_lda %>%
  mutate(
    Comp1 = plsda_df$Comp1,
    Comp2 = plsda_df$Comp2,
    in_tregs_better_region = plsda_df$in_ellipse_plus1,
    in_tregs_worse_region = plsda_df$in_ellipse_minus1,
    region = plsda_df$region,
    treg_outcome = plsda_df$treg_outcome
  )

# Compare parameter values across regions
region_comparison <- df_analysis %>%
  group_by(region) %>%
  summarise(
    n = n(),
    across(all_of(param_names), list(
      mean = mean,
      median = median,
      sd = sd
    ))
  )

# Get just the means for easier viewing
region_means <- df_analysis %>%
  group_by(region) %>%
  summarise(across(all_of(param_names), mean)) %>%
  pivot_longer(cols = -region, names_to = "parameter", values_to = "mean_value") %>%
  pivot_wider(names_from = region, values_from = mean_value)

print(region_means)

# Which parameters differ most between regions?
region_differences <- region_means %>%
  mutate(
    better_vs_worse = abs(`Tregs better only` - `Tregs worse only`),
    better_vs_neutral = abs(`Tregs better only` - `Tregs don't matter`),
    worse_vs_neutral = abs(`Tregs worse only` - `Tregs don't matter`)
  ) %>%
  arrange(desc(better_vs_worse))

print("Parameters most different between Tregs better vs worse regions:")
print(head(region_differences, 10))

# ======== Step 2: Use PLS-DA Loadings to Understand Components
# ===== The PLS-DA loadings tell you which parameters contribute most to each component

# Get PLS-DA loadings
loadings <- plsda_model$loadings$X

# Component 1 interpretation (horizontal axis)
comp1_loadings <- data.frame(
  parameter = rownames(loadings),
  loading = loadings[, 1]
) %>%
  arrange(desc(abs(loading)))

cat("\n========== COMPONENT 1 INTERPRETATION ==========\n")
cat("Component 1 (horizontal axis) is driven by:\n\n")
print(head(comp1_loadings, 10))

# Component 2 interpretation (vertical axis)
comp2_loadings <- data.frame(
  parameter = rownames(loadings),
  loading = loadings[, 2]
) %>%
  arrange(desc(abs(loading)))

cat("\n========== COMPONENT 2 INTERPRETATION ==========\n")
cat("Component 2 (vertical axis) is driven by:\n\n")
print(head(comp2_loadings, 10))

# Visualize loadings
# Biplot showing which parameters drive the separation
p_loadings <- ggplot(comp1_loadings %>% head(10), 
                     aes(x = reorder(parameter, loading), y = loading)) +
  geom_col(aes(fill = loading > 0)) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 10 Parameters Driving Component 1",
       x = "Parameter", y = "Loading") +
  scale_fill_manual(values = c("red", "blue"), guide = "none")

print(p_loadings)

# ==== Step 3: Characterize Each Region
# Statistical tests: which parameters are significantly different?
# Test each parameter between Tregs better vs worse regions
param_tests <- data.frame()

for(param in param_names) {
  # Get values for each region
  better_vals <- df_analysis %>% 
    filter(region == "Tregs better only") %>% 
    pull(!!sym(param))
  
  worse_vals <- df_analysis %>% 
    filter(region == "Tregs worse only") %>% 
    pull(!!sym(param))
  
  neutral_vals <- df_analysis %>% 
    filter(region == "Tregs don't matter") %>% 
    pull(!!sym(param))
  
  # Skip if not enough data
  if(length(better_vals) < 2 || length(worse_vals) < 2) next
  
  # Wilcoxon test
  test_better_vs_worse <- wilcox.test(better_vals, worse_vals)
  test_better_vs_neutral <- wilcox.test(better_vals, neutral_vals)
  test_worse_vs_neutral <- wilcox.test(worse_vals, neutral_vals)
  
  # Effect sizes (Cohen's d approximation)
  d_better_worse <- (mean(better_vals) - mean(worse_vals)) / 
    sqrt((sd(better_vals)^2 + sd(worse_vals)^2) / 2)
  
  param_tests <- rbind(param_tests, data.frame(
    parameter = param,
    mean_better = mean(better_vals),
    mean_worse = mean(worse_vals),
    mean_neutral = mean(neutral_vals),
    p_better_vs_worse = test_better_vs_worse$p.value,
    p_better_vs_neutral = test_better_vs_neutral$p.value,
    p_worse_vs_neutral = test_worse_vs_neutral$p.value,
    cohens_d = d_better_worse
  ))
}

# Adjust for multiple testing
param_tests <- param_tests %>%
  mutate(
    adj_p_better_vs_worse = p.adjust(p_better_vs_worse, method = "BH"),
    adj_p_better_vs_neutral = p.adjust(p_better_vs_neutral, method = "BH"),
    adj_p_worse_vs_neutral = p.adjust(p_worse_vs_neutral, method = "BH")
  ) %>%
  arrange(p_better_vs_worse)

print("Parameters significantly different between regions:")
print(param_tests %>% filter(adj_p_better_vs_worse < 0.05))

# ======== Step 4: Visualize Key Parameters
# Get top discriminating parameters
top_params <- head(param_tests$parameter, 6)

# Boxplots for top parameters
df_plot_params <- df_analysis %>%
  dplyr::filter(region %in% c("Tregs better only", "Tregs worse only", "Tregs don't matter")) %>%
  dplyr::select(region, all_of(top_params)) %>%
  pivot_longer(cols = -region, names_to = "parameter", values_to = "value")

p_params <- ggplot(df_plot_params, aes(x = region, y = value, fill = region)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~parameter, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c(
    "Tregs better only" = "lightblue",
    "Tregs worse only" = "pink",
    "Tregs don't matter" = "gray80"
  )) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Top Parameters Distinguishing Treg Effect Regions",
       x = "", y = "Parameter Value")

print(p_params)

ggsave(
  filename = paste0("./treg_top_params_PLS_DA.png"),
  plot = p_params,
  width = 12,
  height = 8,
  dpi = 300,
  bg='white'
)

dev.off()

# ======= Step 5: Create Parameter Profiles
# Normalize parameters to [0,1] for comparison
df_profiles <- df_analysis %>%
  dplyr::filter(region %in% c("Tregs better only", "Tregs worse only")) %>%
  dplyr::select(region, all_of(param_names)) %>%
  dplyr::group_by(region) %>%
  summarise(across(all_of(param_names), mean)) %>%
  pivot_longer(cols = -region, names_to = "parameter", values_to = "value")

# Normalize each parameter
df_profiles <- df_profiles %>%
  dplyr::group_by(parameter) %>%
  dplyr::mutate(
    value_normalized = (value - min(value)) / (max(value) - min(value))
  ) %>%
  ungroup()

# ======== Step 6: Interpretable Rules

# Create decision rules based on regions
cat("\n========== INTERPRETABLE RULES ==========\n")

# For "Tregs better" region
better_rules <- param_tests %>%
  filter(adj_p_better_vs_worse < 0.05) %>%
  head(5) %>%
  mutate(
    rule = ifelse(mean_better > mean_worse,
                  paste0(parameter, " HIGH (>", round(mean_better, 2), ")"),
                  paste0(parameter, " LOW (<", round(mean_better, 2), ")"))
  )

cat("\nTregs are BENEFICIAL when:\n")
cat(paste(better_rules$rule, collapse = "\nAND "))

# For "Tregs worse" region
worse_rules <- param_tests %>%
  filter(adj_p_better_vs_worse < 0.05) %>%
  head(5) %>%
  mutate(
    rule = ifelse(mean_worse > mean_better,
                  paste0(parameter, " HIGH (>", round(mean_worse, 2), ")"),
                  paste0(parameter, " LOW (<", round(mean_worse, 2), ")"))
  )

cat("\n\nTregs are HARMFUL when:\n")
cat(paste(worse_rules$rule, collapse = "\nAND "))

