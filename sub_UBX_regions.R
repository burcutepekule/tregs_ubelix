dfp    = df_comparisons_plot %>%
  mutate(
    x = abs(cohens_d),
    y = diff_tregs_on_minus_off
  )

# Calculate Euclidean distance from origin (0,0)
dfp = dfp %>%
  mutate(euclidean_dist = sqrt(x^2 + (y/125)^2)) #otherwise y dominates!

# Create region classification
dfp = dfp %>%
  mutate(region = case_when(
    x >= cohens_th & y >= e_th & tregs_better_cohens!=0 ~ "blue_region",  # blue region
    x >= cohens_th & y <= -e_th & tregs_better_cohens!=0 ~ "pink_region", # pink region
    TRUE ~ "other"  # everything else stays gray
  ))

saveRDS(dfp,'dfp_merged.rds')

# Create separate color assignment for blue and pink regions
# Use simple direct mapping based on euclidean distance

# For blue region
blue_points = dfp %>% filter(region == "blue_region")
if(nrow(blue_points) > 0) {
  blue_palette = colorRampPalette(c("#637AB9", "#4FB7B3", "#016B61", "#2F5755"))(nrow(blue_points))
  # blue_palette = colorRampPalette(c("#313647", "#2F5755", "#016B61", "#4FB7B3", "#637AB9", "#3B38A0", "#1A2A80"))(100)
  blue_dist_breaks = seq(min(blue_points$euclidean_dist), max(blue_points$euclidean_dist), length.out = nrow(blue_points)+1)
  blue_points$color_idx = cut(blue_points$euclidean_dist, breaks = blue_dist_breaks, labels = FALSE, include.lowest = TRUE)
  blue_points$point_color = blue_palette[blue_points$color_idx]
}

# For pink region  
pink_points = dfp %>% filter(region == "pink_region")
if(nrow(pink_points) > 0) {
  red_palette = colorRampPalette(c("#FF9B17", "#FA812F", "#E83F25", "#BF092F"))(nrow(pink_points))
  pink_dist_breaks = seq(min(pink_points$euclidean_dist), max(pink_points$euclidean_dist), length.out = nrow(pink_points)+1)
  pink_points$color_idx = cut(pink_points$euclidean_dist, breaks = pink_dist_breaks, labels = FALSE, include.lowest = TRUE)
  pink_points$point_color = red_palette[pink_points$color_idx]
}

# For other region
other_points = dfp %>% filter(region == "other")
other_points$point_color = "gray70"

# Combine all
dfp_final = bind_rows(blue_points, pink_points, other_points)

# Create plot
p = ggplot(dfp_final, aes(x = x, y = y, shape = injury_type, color = point_color)) +
  annotate("rect",
           xmin = cohens_th, xmax = Inf,
           ymin = -Inf, ymax = -1*e_th,
           fill = "pink", alpha = 0.3) +
  annotate("rect",
           xmin = cohens_th, xmax = Inf,
           ymin = e_th, ymax = Inf,
           fill = "lightblue", alpha = 0.3) +
  annotate("rect",
           xmin = -Inf, xmax = cohens_th,
           ymin = -Inf, ymax = Inf,
           fill = "gray80", alpha = 0.15) +
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = -1*e_th, ymax = e_th,
           fill = "gray80", alpha = 0.15) +
  geom_point(size = 3) +
  scale_shape_manual(name = "Injury Type", values = c("sterile" = 16, "pathogenic" = 2)) +
  scale_color_identity(guide = "none") +  # Use colors directly
  geom_vline(xintercept = cohens_th, linetype = "dashed") +
  geom_hline(yintercept = e_th, linetype = "dashed") +
  geom_hline(yintercept = -1*e_th, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "solid", col='black', size = 0.7) +
  theme_minimal() +
  labs(x = paste0("jensen-shannon distance (threshold = ", cohens_th, ")"),
       y = paste0("Epithelial score for tregs_on - tregs_off (threshold = ±", round(e_th), " for max score of 125.)"),
       title = paste0('Num of parameter sets: ',length(unique(dfp$param_set_id))),
       shape = "Injury Type") + 
  scale_y_continuous(
    breaks = sort(c(seq(round(min(dfp$diff_tregs_on_minus_off))-5,
                        round(max(dfp$diff_tregs_on_minus_off))+5, by=10), 0))
  ) +
  scale_x_continuous(
    breaks = sort(c(seq(0, 1, by=0.1), cohens_th))
  ) +
  guides(
    shape = guide_legend(order = 1)
  )

