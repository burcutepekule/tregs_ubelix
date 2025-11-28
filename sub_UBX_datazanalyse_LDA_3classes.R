library(MASS)
library(dplyr)
library(ggplot2)

# LDA model (using same data as PLS-DA)
lda_model = lda(
  x = df_lda %>% dplyr::select(all_of(param_names)),
  grouping = df_lda$tregs_better_cohens
)

# # Use regularized LDA (which handles collinearity)
# library(klaR)
# lda_model <- rda(tregs_better_cohens ~ ., 
#                   data = df_lda %>% dplyr::select(all_of(param_names), tregs_better_cohens),
#                   gamma = 0.5, lambda = 0.5)


# Get LDA projections
lda_pred = predict(lda_model)
lda_df = data.frame(
  LD1 = lda_pred$x[, 1],
  LD2 = if(ncol(lda_pred$x) > 1) lda_pred$x[, 2] else 0,  # Handle case with only 1 LD
  treg_outcome = as.factor(df_lda$tregs_better_cohens)
)

# ============ DUAL ELLIPSE ANALYSIS ============

# Function to analyze ellipse for a specific class (same as before)
analyze_class_ellipse = function(lda_df, target_class, conf_levels) {
  
  # Get points for target class
  class_points = lda_df %>% filter(treg_outcome == target_class)
  
  # Fit ellipse to class distribution
  class_center = colMeans(class_points[, c("LD1", "LD2")])
  class_cov = cov(class_points[, c("LD1", "LD2")])
  
  # Calculate Mahalanobis distance for all points
  mahal_dist = mahalanobis(
    lda_df[, c("LD1", "LD2")],
    center = class_center,
    cov = class_cov
  )
  
  ellipse_results = data.frame()
  
  for(conf_level in conf_levels) {
    chi2_threshold = qchisq(conf_level, df = 2)
    in_ellipse = mahal_dist < chi2_threshold
    
    # Calculate metrics
    tp = sum(lda_df$treg_outcome == target_class & in_ellipse)
    fp = sum(lda_df$treg_outcome != target_class & in_ellipse)
    fn = sum(lda_df$treg_outcome == target_class & !in_ellipse)
    
    total_target = sum(lda_df$treg_outcome == target_class)
    total_other = sum(lda_df$treg_outcome != target_class)
    
    recall = tp / total_target
    precision = tp / (tp + fp)
    pct_other_in = fp / total_other
    enrichment = recall / pct_other_in
    
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
conf_levels = c(0.50, 0.75, 0.90, 0.95, 0.99)

ellipse_plus1_lda = analyze_class_ellipse(lda_df, "1", conf_levels)
ellipse_minus1_lda = analyze_class_ellipse(lda_df, "-1", conf_levels)

# ============ VISUALIZATIONS ============

# Choose confidence levels for visualization
# level_plus1  = 0.75
# level_minus1 = 0.75

# Add ellipse membership to dataframe
lda_df$in_ellipse_plus1  = ellipse_plus1_lda$mahal_dist < qchisq(level_plus1, df = 2)
lda_df$in_ellipse_minus1 = ellipse_minus1_lda$mahal_dist < qchisq(level_minus1, df = 2)

# Create a combined category
lda_df = lda_df %>%
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
  
  radius = sqrt(qchisq(level, df = 2))
  ellipse = circle %*% chol(cov_matrix) * radius
  ellipse = sweep(ellipse, 2, center, "+")
  
  return(data.frame(x = ellipse[, 1], y = ellipse[, 2]))
}

# Generate ellipse coordinates
ellipse_plus1_coords_lda = get_ellipse_coords(
  ellipse_plus1_lda$center, 
  ellipse_plus1_lda$cov, 
  level_plus1
)

ellipse_minus1_coords_lda = get_ellipse_coords(
  ellipse_minus1_lda$center, 
  ellipse_minus1_lda$cov, 
  level_minus1
)

# ============ COMPARISON SUMMARY ============
cat("\n========== LDA vs PLS-DA COMPARISON ==========\n")
cat("\nLDA Results:\n")
cat(sprintf("  Tregs better enrichment: %.1f-fold\n", 
            ellipse_plus1_lda$results$enrichment[ellipse_plus1_lda$results$conf_level == level_plus1]))
cat(sprintf("  Tregs worse enrichment: %.1f-fold\n", 
            ellipse_minus1_lda$results$enrichment[ellipse_minus1_lda$results$conf_level == level_minus1]))

# Print LDA loadings (coefficients of linear discriminants)
cat("\n========== LDA LOADINGS (LD1) ==========\n")
lda_loadings_ld1 = data.frame(
  parameter = rownames(lda_model$scaling),
  loading = lda_model$scaling[, 1]
) %>%
  arrange(desc(abs(loading)))

print(head(lda_loadings_ld1, 10))

if(ncol(lda_model$scaling) > 1) {
  cat("\n========== LDA LOADINGS (LD2) ==========\n")
  lda_loadings_ld2 = data.frame(
    parameter = rownames(lda_model$scaling),
    loading = lda_model$scaling[, 2]
  ) %>%
    arrange(desc(abs(loading)))
  
  print(head(lda_loadings_ld2, 10))
}


# Create a scaling factor for the arrows so they're visible
arrow_scale = 1  # Adjust this to make arrows longer/shorter

# Get the loadings and scale them
lda_arrows = data.frame(
  parameter = rownames(lda_model$scaling),
  LD1 = lda_model$scaling[, 1] * arrow_scale,
  LD2 = if(ncol(lda_model$scaling) > 1) lda_model$scaling[, 2] * arrow_scale else 0
)

# Select top N most important parameters (by total contribution)
n_arrows = 12 #23 to select all 
lda_arrows = lda_arrows %>%
  mutate(total_loading = sqrt(LD1^2 + LD2^2)) %>%
  arrange(desc(total_loading)) %>%
  head(n_arrows)

# Plot with arrows
p2_lda = ggplot(lda_df, aes(x = LD1, y = LD2)) +
  geom_point(data = lda_df %>% filter(treg_outcome == "0"),
             aes(color = region), alpha = 0.8, size = 2) +
  scale_color_manual(values = c(
    "Tregs don't matter" = "gray90",
    "Tregs better only" = "lightblue",
    "Tregs worse only" = "pink",
    "Both" = "orange"
  )) +
  geom_path(data = ellipse_plus1_coords_lda, aes(x = x, y = y),
            color = "blue", linewidth = 1.5, inherit.aes = FALSE) +
  geom_path(data = ellipse_minus1_coords_lda, aes(x = x, y = y),
            color = "red", linewidth = 1.5, inherit.aes = FALSE) +
  geom_point(data = lda_df %>% filter(treg_outcome == "1"),
             color = "darkblue", alpha = 0.8, size = 3) +
  geom_point(data = lda_df %>% filter(treg_outcome == "-1"),
             color = "darkred", alpha = 0.8, size = 3) +
  # Add parameter vectors
  geom_segment(data = lda_arrows,
               aes(x = 0, y = 0, xend = LD1, yend = LD2),
               arrow = arrow(length = unit(0.3, "cm")),
               color = "black", linewidth = 0.4,
               inherit.aes = FALSE) +
  geom_text_repel(data = lda_arrows,
                  aes(x = LD1, y = LD2, label = parameter),
                  color = "black", size = 3, fontface = "bold",
                  inherit.aes = FALSE) +
  theme_minimal() +
  labs(title = paste0("LDA of parameter sets based on Treg help: conf. : (",level_plus1,", ",level_minus1,")"),
       subtitle = sprintf("Enrichment Tregs better: %.1f-fold | Enrichment Tregs worse: %.1f-fold", 
                          ellipse_plus1_lda$results$enrichment[ellipse_plus1_lda$results$conf_level == level_plus1],
                          ellipse_minus1_lda$results$enrichment[ellipse_minus1_lda$results$conf_level == level_minus1]))

print(p2_lda)

ggsave(
  filename = paste0("./treg_regions_LDA_",inj_type,"_w_arrows.png"),
  plot = p2_lda,
  width = 9,
  height = 6,
  dpi = 300,
  bg='white'
)

dev.off()