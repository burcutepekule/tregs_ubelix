# ========= Step 1: Extract Parameter Values for Each ***cohens category*** (region is defined over that)
# Add region membership and original parameters to analysis dataframe
df_analysis = df_lda %>%
  mutate(
    Comp1 = plsda_df$Comp1,
    Comp2 = plsda_df$Comp2,
    in_tregs_better_region = (tregs_better_cohens==1),
    in_tregs_worse_region = (tregs_better_cohens==-1),
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
  
  worse_vals = df_analysis %>% 
    filter(region == "Tregs worse only") %>%
    pull(!!sym(param))
  
  # Skip if not enough data
  if(length(better_vals)>2 & length(neutral_vals)>2 & length(worse_vals)<2){
    # Wilcoxon test
    test_better_vs_neutral = wilcox.test(better_vals, neutral_vals)
    # Adjust for multiple testing
    param_tests = rbind(param_tests, data.frame(
      p_better_vs_neutral = test_better_vs_neutral$p.value,
      adj_p_better_vs_neutral = p.adjust(test_better_vs_neutral$p.value, method = "BH")
    ))
  }else if(length(worse_vals)>2 & length(neutral_vals)>2 & length(better_vals)<2){
    # Wilcoxon test
    p_worse_vs_neutral = test_worse_vs_neutral$p.value
    # Adjust for multiple testing
    param_tests = rbind(param_tests, data.frame(
      p_worse_vs_neutral = test_worse_vs_neutral$p.value,
      adj_p_worse_vs_neutral = p.adjust(test_worse_vs_neutral$p.value, method = "BH")
    ))
  }else if(length(better_vals)>2 & length(worse_vals)>2 & length(neutral_vals)<2){
    # Wilcoxon test
    test_better_vs_worse = wilcox.test(better_vals, worse_vals)
    # Adjust for multiple testing
    param_tests = rbind(param_tests, data.frame(
      p_better_vs_worse = test_better_vs_worse$p.value,
      adj_p_better_vs_worse = p.adjust(test_better_vs_worse$p.value, method = "BH")
    )) 
  }else{
    # Wilcoxon test
    test_better_vs_worse = wilcox.test(better_vals, worse_vals)
    test_better_vs_neutral = wilcox.test(better_vals, neutral_vals)
    test_worse_vs_neutral = wilcox.test(worse_vals, neutral_vals)
    # Adjust for multiple testing
    param_tests = rbind(param_tests, data.frame(
      p_better_vs_worse = test_better_vs_worse$p.value,
      p_better_vs_neutral = test_better_vs_neutral$p.value,
      p_worse_vs_neutral = test_worse_vs_neutral$p.value,
      adj_p_better_vs_worse = p.adjust(test_better_vs_worse$p.value, method = "BH"),
      adj_p_better_vs_neutral = p.adjust(test_better_vs_neutral$p.value, method = "BH"),
      adj_p_worse_vs_neutral = p.adjust(test_worse_vs_neutral$p.value, method = "BH")
    )) 
  }
}

# ======== Step 4: Visualize Key Parameters

# Boxplots for top parameters
df_plot_params = df_analysis %>%
  dplyr::select(region, all_of(param_names)) %>%
  pivot_longer(cols = -region, names_to = "parameter", values_to = "value") 

# Update y_positions
y_positions = df_plot_params %>%
  group_by(parameter) %>%
  summarise(y_pos = max(value, na.rm = TRUE) * 1.1, .groups = "drop")

regions_plot  = unique(df_plot_params$region)
regions_combo = combn(regions_plot, 2, simplify = FALSE)


if(violin_on==1){
  p_params = ggplot(df_plot_params, aes(x = region, y = value, fill = region)) +
    geom_violin(alpha = 0.2, trim = TRUE) +
    geom_boxplot(width = 0.2, alpha = 0.2) +
    # facet_wrap(~parameter_labeled, scales = "free_y", ncol = 4) + # if you want it ordered according to p val
    facet_wrap(~parameter, scales = "free_y", ncol = 4) + # ordered alphabetically (easier to compare among inj types.)
    scale_fill_manual(values = c(
      "Tregs better only" = "lightblue",
      "Tregs don't matter" = "gray80",
      "Tregs worse only" = "pink"
    )) +
    scale_y_continuous(expand = expansion(mult = c(0.00, 0.10))) +  # Add this line
    # Add significance brackets
    geom_signif(
      comparisons = regions_combo,
      test = "t.test",
      test.args = list(var.equal = FALSE),   # Welch t-test
      map_signif_level = TRUE,
      textsize = 5,
      step_increase = 0.1,
      tip_length = 0.02
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 12)
    ) +
    labs(title = "Top Parameters Distinguishing Treg Effect Regions",
         x = "", y = "Parameter Value")
  
}else{
  
  p_params = ggplot(df_plot_params, aes(x = region, y = value, fill = region)) +
    geom_violin(alpha = 0.2, trim = TRUE) +
    geom_boxplot(width = 0.2, alpha = 0.2) +
    # facet_wrap(~parameter_labeled, scales = "free_y", ncol = 4) + # if you want it ordered according to p val
    facet_wrap(~parameter, scales = "free_y", ncol = 4) + # ordered alphabetically (easier to compare among inj types.)
    scale_fill_manual(values = c(
      "Tregs better only" = "lightblue",
      "Tregs don't matter" = "gray80",
      "Tregs worse only" = "pink"
    )) +
    scale_y_continuous(expand = expansion(mult = c(0.00, 0.10))) +  # Add this line
    # Add significance brackets
    geom_signif(
      comparisons = regions_combo,
      test = "t.test",
      test.args = list(var.equal = FALSE),   # Welch t-test
      map_signif_level = TRUE,
      textsize = 5,
      step_increase = 0.1,
      tip_length = 0.02
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 12)
    ) +
    labs(title = "Top Parameters Distinguishing Treg Effect Regions",
         x = "", y = "Parameter Value")
  
}

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
