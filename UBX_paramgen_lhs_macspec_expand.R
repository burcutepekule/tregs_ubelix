rm(list=ls())
library(dplyr)
library(tidyr)
library(lhs)
library(readr)

setwd('/Users/burcutepekule/Dropbox/tregs_ubelix/')

# == The csv file below and dfp dataframe was saved by  
# == ~/Dropbox/tregs_ubelix/UBX_datazanalyse_jsd_regions_merged_dataset_nolabel.R
# == and ~/Dropbox/tregs_ubelix/sub_UBX_regions.R respectively
df_params = read_csv('/Users/burcutepekule/Dropbox/tregs_ubelix/lhs_parameters_ubelix_merged.csv', show_col_types = FALSE)
dfp       = readRDS('dfp_merged.rds')
ids_tregs_better = unique(dfp %>% dplyr::filter(region=='blue_region') %>% pull(param_set_id))
ids_tregs_worse  = unique(dfp %>% dplyr::filter(region=='pink_region') %>% pull(param_set_id))
ids_tregs_dm     = unique(dfp %>% dplyr::filter(region=='other') %>% pull(param_set_id))

# === add the new parameters
df_params$mac_discrimination_efficiency = 1 # compare with the perfect scenario!
N = 9 # 
df_params_expanded = df_params %>%
  slice(rep(1:n(), each = N)) %>%
  dplyr::group_by(param_set_id) %>%
  dplyr::mutate(
    mac_rat_com_pat_threshold = seq(0.5, 1, length.out = N),
    param_set_id_new = paste0(param_set_id, "_", row_number())
  ) %>%
  ungroup()

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
  "mac_discrimination_efficiency",
  "mac_rat_com_pat_threshold"
)
# === add the new parameters

df_params_tregs_better = df_params_expanded %>% dplyr::filter(param_set_id %in% ids_tregs_better)
df_params_tregs_worse  = df_params_expanded %>% dplyr::filter(param_set_id %in% ids_tregs_worse)
df_params_tregs_dm     = df_params_expanded %>% dplyr::filter(param_set_id %in% ids_tregs_dm)

# Reorder columns to put param_set_id first
samples_better = df_params_tregs_better[c('param_set_id_new', param_names)]
samples_worse  = df_params_tregs_worse[c('param_set_id_new', param_names)]
samples_dm     = df_params_tregs_dm[c('param_set_id_new', param_names)]

# Export
write.csv(samples_better, "./lhs_parameters_macspec_expanded_tregs_better.csv", row.names = FALSE)
write.csv(samples_worse, "./lhs_parameters_macspec_expanded_tregs_worse.csv", row.names = FALSE)
write.csv(samples_dm, "./lhs_parameters_macspec_expanded_tregs_dm.csv", row.names = FALSE)

