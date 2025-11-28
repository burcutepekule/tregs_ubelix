# ========= Step 1: Extract Parameter Values for Each ***cohens category*** (region is defined over that)
# Add region membership and original parameters to analysis dataframe
df_analysis = df_lda %>%
  mutate(
    Comp1 = plsda_df$Comp1,
    Comp2 = plsda_df$Comp2,
    in_tregs_better_region = (tregs_better_cohens==1),
    in_tregs_worse_region =  (tregs_better_cohens==-1),
    region = ifelse(tregs_better_cohens==-1, "Tregs worse only", 
                    ifelse(tregs_better_cohens==1, "Tregs better only","Tregs don't matter")),
    treg_outcome = tregs_better_cohens
  )

# Compare parameter values across regions
region_comparison = df_analysis %>%
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
region_means = df_analysis %>%
  group_by(region) %>%
  summarise(across(all_of(param_names), mean)) %>%
  pivot_longer(cols = -region, names_to = "parameter", values_to = "mean_value") %>%
  pivot_wider(names_from = region, values_from = mean_value)

print(region_means)

# Which parameters differ most between regions?
region_differences = region_means %>%
  mutate(
    better_vs_neutral = abs(`Tregs better only` - `Tregs don't matter`)/(`Tregs better only`+`Tregs don't matter`)
  ) %>%
  arrange(desc(better_vs_neutral))

print("Parameters most different between Tregs better vs worse regions:")
print(head(region_differences, 10))

# ======== Step 2: Use PLS-DA Loadings to Understand Components
# ===== The PLS-DA loadings tell you which parameters contribute most to each component

# Get PLS-DA loadings
loadings = plsda_model$loadings$X

# Component 1 interpretation (horizontal axis)
comp1_loadings = data.frame(
  parameter = rownames(loadings),
  loading = loadings[, 1]
) %>%
  arrange(desc(abs(loading)))

# cat("\n========== COMPONENT 1 INTERPRETATION ==========\n")
# cat("Component 1 (horizontal axis) is driven by:\n\n")
# print(head(comp1_loadings, 10))

# Component 2 interpretation (vertical axis)
comp2_loadings = data.frame(
  parameter = rownames(loadings),
  loading = loadings[, 2]
) %>%
  arrange(desc(abs(loading)))

# cat("\n========== COMPONENT 2 INTERPRETATION ==========\n")
# cat("Component 2 (vertical axis) is driven by:\n\n")
# print(head(comp2_loadings, 10))

# Visualize loadings
# Biplot showing which parameters drive the separation
p_loadings = ggplot(comp1_loadings %>% head(10), 
                     aes(x = reorder(parameter, loading), y = loading)) +
  geom_col(aes(fill = loading > 0)) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 10 Parameters Driving Component 1",
       x = "Parameter", y = "Loading") +
  scale_fill_manual(values = c("red", "blue"), guide = "none")

# print(p_loadings)

# ==== Step 3: Characterize Each Region
# Statistical tests: which parameters are significantly different?
# Test each parameter between Tregs better vs worse regions
param_tests = data.frame()

for(param in param_names) {
  # Get values for each region
  # Not region but the real classification
  better_vals = df_analysis %>% 
    filter(region == "Tregs better only") %>%
    pull(!!sym(param))
  
  neutral_vals = df_analysis %>% 
    filter(region == "Tregs don't matter") %>%
    pull(!!sym(param))
  
  # Skip if not enough data
  if(length(better_vals) < 2) next
  
  # Wilcoxon test
  test_better_vs_neutral = wilcox.test(better_vals, neutral_vals)

  param_tests = rbind(param_tests, data.frame(
    parameter = param,
    mean_better = mean(better_vals),
    mean_neutral = mean(neutral_vals),
    p_better_vs_neutral = test_better_vs_neutral$p.value
  ))
}

# Adjust for multiple testing
param_tests = param_tests %>%
  mutate(
    adj_p_better_vs_neutral = p.adjust(p_better_vs_neutral, method = "BH")
  ) %>%
  arrange(p_better_vs_neutral)

print("Parameters significantly different between regions:")
print(param_tests %>% filter(adj_p_better_vs_neutral < 0.05))

# ======== Step 4: Visualize Key Parameters
# Get top discriminating parameters
# top_params = head(param_tests$parameter, 10)
top_params = param_tests$parameter

# Get p-values for the top parameters
p_values_for_plot = param_tests %>%
  filter(parameter %in% top_params) %>%
  dplyr::select(parameter, adj_p_better_vs_neutral) %>%
  mutate(
    p_label_2 = case_when(
      adj_p_better_vs_neutral < 0.001 ~ "***",
      adj_p_better_vs_neutral < 0.01 ~ "**",
      adj_p_better_vs_neutral < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    p_text = sprintf("p_bvn = %.2e %s", adj_p_better_vs_neutral, p_label_2)
  )

# Boxplots for top parameters
df_plot_params = df_analysis %>%
  dplyr::filter(region %in% c("Tregs better only", "Tregs don't matter")) %>%
  dplyr::select(region, all_of(top_params)) %>%
  pivot_longer(cols = -region, names_to = "parameter", values_to = "value") %>%
  mutate(
    param_order = match(parameter, top_params),
    parameter_labeled = paste0(sprintf("%02d", param_order), ". ", parameter),
    parameter_labeled = factor(parameter_labeled, 
                               levels = paste0(sprintf("%02d", 1:length(top_params)), ". ", top_params))
  )

# Update y_positions
y_positions = df_plot_params %>%
  group_by(parameter, parameter_labeled) %>%
  summarise(y_pos = max(value, na.rm = TRUE) * 1.1, .groups = "drop") %>%
  left_join(p_values_for_plot, by = "parameter")

p_params = ggplot(df_plot_params, aes(x = region, y = value, fill = region)) +
  geom_violin(alpha = 0.2, trim = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
  # facet_wrap(~parameter_labeled, scales = "free_y", ncol = 4) + # if you want it ordered according to p val
  facet_wrap(~parameter, scales = "free_y", ncol = 4) + # ordered alphabetically (easier to compare among inj types.)
  scale_fill_manual(values = c(
    "Tregs better only" = "lightblue",
    "Tregs don't matter" = "gray80"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.00, 0.10))) +  # Add this line
  # Add significance brackets
  # geom_signif(
  #   comparisons = list(
  #     c("Tregs better only", "Tregs don't matter"),
  #     textsize = 5,
  #     step_increase = 0.2,
  #     tip_length = 0.02
  #   ),
  #   map_signif_level = TRUE,
  #   textsize = 5,
  #   step_increase = 0.1,  # Vertical spacing between brackets
  #   tip_length = 0.02     # Length of bracket tips
  # ) +
  geom_signif(
    comparisons = list(c("Tregs better only", "Tregs don't matter")),
    test = "t.test",
    test.args = list(var.equal = FALSE),   # Welch t-test
    map_signif_level = TRUE,
    textsize = 5,
    step_increase = 0.1,
    tip_length = 0.02
  )+
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 12)
  ) +
  labs(title = "Top Parameters Distinguishing Treg Effect Regions",
       x = "", y = "Parameter Value")

print(p_params)

ggsave(
  filename = paste0("./treg_top_params_PLS_DA_",inj_type,"_classification.png"),
  plot = p_params,
  width = 18,
  height = 20,
  dpi = 300,
  bg='white'
)

dev.off()


df_plot_params_pick = df_plot_params %>% dplyr::filter(parameter=='treg_discrimination_efficiency')
table(df_plot_params_pick$region)

df_plot_params_pick_dm = df_plot_params_pick %>% dplyr::filter(region=="Tregs don't matter") %>% dplyr::pull(value)
df_plot_params_pick_bt = df_plot_params_pick %>% dplyr::filter(region=="Tregs better only") %>% dplyr::pull(value)
t.test(df_plot_params_pick_dm, df_plot_params_pick_bt)
