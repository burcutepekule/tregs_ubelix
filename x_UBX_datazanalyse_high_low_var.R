rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(readr)  # For read_csv
library(stringr)
library(zoo)
library(scales)
library(ggrepel)

source("./MISC/PLOT_FUNCTIONS.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

df_params       = read_csv('/Users/burcutepekule/Dropbox/tregs_ubelix/lhs_parameters_ubelix.csv', show_col_types = FALSE)
df_results_keep = readRDS('/Users/burcutepekule/Dropbox/tregs_ubelix/data_cpp_read_cohens.rds')
length(unique(df_results_keep$param_set_id))

# --- filter for complete # of reps 
reps_df       = as.data.frame(table(df_results_keep$param_set_id))
keep_param_id = reps_df %>% dplyr::filter(Freq==600) %>% dplyr::pull(Var1) # 600 = 100 reps per scenario, 6 scenarios
df_results    = df_results_keep %>% filter(param_set_id %in% keep_param_id)
unique(table(df_results$param_set_id)) #66, sanity check
length(unique(df_results$param_set_id))

#----- filter based on ss_start, it cannot be too large otherwise not much to compare!
ss_start_threshold = 4500
param_id_all_below = df_results %>%
  dplyr::group_by(param_set_id) %>%
  dplyr::summarise(all_below = all(ss_start < ss_start_threshold), .groups = "drop") %>%
  dplyr::filter(all_below) %>%
  dplyr::pull(param_set_id)
length(param_id_all_below)/length(unique(df_results$param_set_id)) # >99%!
df_results = df_results %>% dplyr::filter(param_set_id %in% param_id_all_below)
max(df_results$ss_start)<ss_start_threshold # TRUE, sanity check
unique(table(df_results$param_set_id)) #600, sanity check
length(unique(df_results$param_set_id))

df_averaged = df_results %>%
  group_by(param_set_id, injury_type, tregs_on, tregs_rnd) %>%
  summarise(mean_mean_score = mean(mean_score), 
            sd_mean_score = sd(mean_score),
            .groups = 'drop')

df_comparisons = df_averaged %>%
  # Create a key combining the conditions
  mutate(condition = case_when(
    tregs_on == 0 & tregs_rnd == 0 ~ "tregs_off",
    tregs_on == 1 & tregs_rnd == 0 ~ "tregs_on_no_rnd",
    tregs_on == 1 & tregs_rnd == 1 ~ "tregs_on_rnd"
  )) %>%
  select(param_set_id, injury_type, condition, mean_mean_score) %>%
  pivot_wider(names_from = condition, values_from = mean_mean_score) %>%
  mutate(
    diff_tregs_on_minus_off = tregs_on_no_rnd - tregs_off,
    diff_notrnd_minus_rnd   = tregs_on_no_rnd - tregs_on_rnd
  )
#----- filter based on cohen's d since you are looking for differences!
tol_in              = 25*0.05
cohens_d_th         = 0.5


hist(df_comparisons$diff_tregs_on_minus_off, 30)
min(df_comparisons$diff_tregs_on_minus_off)
hist(df_comparisons$diff_notrnd_minus_rnd, 30)
min(df_comparisons$diff_notrnd_minus_rnd)

df_comparisons_sense    = merge(df_comparisons, distinct(df_results[c('param_set_id','d_10','d_21','d_43','d_54')]), by='param_set_id')
df_averaged             = df_averaged %>% dplyr::mutate(high_var = ifelse(sd_mean_score>=1,1,0))
df_averaged             = df_averaged %>% inner_join(df_params, by='param_set_id')

# ---- HERE: Pick your battle - which case are you checking?
df_plot = df_averaged %>% dplyr::filter(injury_type=='sterile' & tregs_on==1 & tregs_rnd==0) %>%
  dplyr::filter(high_var %in% c(0, 1)) %>%
  dplyr::mutate(high_var = factor(high_var, labels = c("Low Variance", "High Variance")))

df_plot = df_plot %>% dplyr::mutate(activity_engulf_M1_M2_diff = activity_engulf_M1_baseline-activity_engulf_M2_baseline)
# vector of parameter names
param_names = c(
  "th_ROS_microbe",
  "th_ROS_epith_recover",
  "epith_recovery_chance",
  "rat_com_pat_threshold",
  "diffusion_speed_DAMPs",
  "diffusion_speed_SAMPs",
  "diffusion_speed_ROS",
  "add_ROS",
  "add_DAMPs",
  "add_SAMPs",
  "ros_decay",
  "DAMPs_decay",
  "SAMPs_decay",
  "activation_threshold_DAMPs",
  "activation_threshold_SAMPs",
  "activity_engulf_M0_baseline",
  "activity_engulf_M1_baseline",
  "activity_engulf_M2_baseline",
  "activity_ROS_M1_baseline",
  "rate_leak_commensal_injury",
  "rate_leak_pathogen_injury",
  "rate_leak_commensal_baseline",
  "active_age_limit",
  "treg_discrimination_efficiency",
  "activity_engulf_M1_M2_diff"
)

# loop over parameters
dir_name = '/Users/burcutepekule/Dropbox/tregs_ubelix/high_low_var_histograms'
dir.create(dir_name, showWarnings = FALSE)
for (param in param_names) {
  p = ggplot(df_plot, aes(x = .data[[param]], fill = high_var)) +
    geom_density(alpha = 0.5, bw = 0.05) +
    scale_fill_manual(values = c("Low Variance" = "blue", "High Variance" = "red")) +
    labs(x = param, title = paste("Density of", param)) +
    theme_minimal()
  ggsave(
    filename = paste0(dir_name, "/sterile_high_var_", param, "_density.png"),
    plot = p,
    width = 9,
    height = 6,
    dpi = 300
  )
}

