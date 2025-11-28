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

if(any(is.na(ellipse_minus1$center)) & any(!is.na(ellipse_plus1$center))){ #only plus 1
  # ============ VISUALIZATIONS ============
  
  # Add ellipse membership to dataframe
  plsda_df$in_ellipse_plus1  = ellipse_plus1$mahal_dist < qchisq(level_plus1, df = 2)

  # Create a combined category
  plsda_df = plsda_df %>%
    mutate(region = case_when(
      in_ellipse_plus1 ~ "Tregs better only",
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
      "Tregs better only" = "lightblue"
    )) +
    # Manual ellipses using exact calculations
    geom_path(data = ellipse_plus1_coords, aes(x = x, y = y),
              color = "blue", linewidth = 1.5, inherit.aes = FALSE) +
    geom_point(data = plsda_df %>% filter(treg_outcome == "1"),
               color = "darkblue", alpha = 0.8, size = 3) +
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
    labs(title = paste0("PLS DA of parameter sets based on Treg help: conf. : (",level_plus1,")"),
         subtitle = sprintf("Enrichment Tregs better: %.1f-fold", 
                            ellipse_plus1$results$enrichment[ellipse_plus1$results$conf_level == level_plus1]))
  
  
}else if(any(!is.na(ellipse_minus1$center)) & any(is.na(ellipse_plus1$center))){ #only minus 1
  # ============ VISUALIZATIONS ============
  
  # Add ellipse membership to dataframe
  plsda_df$in_ellipse_minus1 = ellipse_minus1$mahal_dist < qchisq(level_minus1, df = 2)
  
  # Create a combined category
  plsda_df = plsda_df %>%
    mutate(region = case_when(
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
      "Tregs worse only" = "pink"
    )) +
    # Manual ellipses using exact calculations
    geom_path(data = ellipse_minus1_coords, aes(x = x, y = y),
              color = "red", linewidth = 1.5, inherit.aes = FALSE) +
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
    labs(title = paste0("PLS DA of parameter sets based on Treg help: conf. : (",level_minus1,")"),
         subtitle = sprintf("Enrichment Tregs worse: %.1f-fold", 
                            ellipse_minus1$results$enrichment[ellipse_minus1$results$conf_level == level_minus1]))
  
  
}else{ #both exist
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
  
  
}

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
