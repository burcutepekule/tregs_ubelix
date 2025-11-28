rm(list=ls())
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

setwd('/Users/burcutepekule/Dropbox/tregs_ubelix')
source("./MISC/PLOT_FUNCTIONS.R")
source("./MISC/DATA_READ_FUNCTIONS.R")


colnames_insert_fixed = c("t","sterile","tregs_on","randomize_tregs","param_set_id","rep_id","epithelial_score","time_ss",
                          'epithelial_healthy','epithelial_inj_1','epithelial_inj_2','epithelial_inj_3','epithelial_inj_4',
                          'epithelial_inj_5','phagocyte_M0','phagocyte_M1_L_0','phagocyte_M1_L_1','phagocyte_M1_L_2','phagocyte_M1_L_3',
                          'phagocyte_M1_L_4','phagocyte_M1_L_5','phagocyte_M2_L_0','phagocyte_M2_L_1','phagocyte_M2_L_2',
                          'phagocyte_M2_L_3','phagocyte_M2_L_4','phagocyte_M2_L_5','commensal','pathogen','treg_resting',
                          'treg_active','C_ROS','C_M0','C_M1','C_M2','P_ROS','P_M0','P_M1','P_M2')

path         = "/Users/burcutepekule/Desktop/mass_sim_results_R_cpp/"

param_id     = 59701
param_id     = 36802
param_id     = 29510
param_id     = 67500

results_0_0_0 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', param_id, '_sterile_0_tregs_0_trnd_0.rds'))
results_0_1_0 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', param_id, '_sterile_0_tregs_1_trnd_0.rds'))
results_0_1_1 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', param_id, '_sterile_0_tregs_1_trnd_1.rds'))

results_1_0_0 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', param_id, '_sterile_1_tregs_0_trnd_0.rds'))
results_1_1_0 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', param_id, '_sterile_1_tregs_1_trnd_0.rds'))
results_1_1_1 = readRDS(paste0(path, 'longitudinal_df_param_set_id_', param_id, '_sterile_1_tregs_1_trnd_1.rds'))

results = rbind(
  results_0_0_0,
  results_0_1_0,
  results_0_1_1,
  results_1_0_0,
  results_1_1_0,
  results_1_1_1
)

colnames(results) = colnames_insert_fixed # this is fixing the mistake in UBX_datagen_cpp.R -> 7:37 should have been 9:39

colnames(results_0_1_0) = colnames_insert_fixed
colnames(results_0_0_0) = colnames_insert_fixed

results_0_1_0 = results_0_1_0 %>%
  dplyr::group_by(rep_id) %>%
  dplyr::filter(t >= time_ss) %>%
  dplyr::ungroup()

results_0_0_0 = results_0_0_0 %>%
  dplyr::group_by(rep_id) %>%
  dplyr::filter(t >= time_ss) %>%
  dplyr::ungroup()


# Combine the data with a grouping variable
combined_data = data.frame(
  epithelial_score = c(results_0_1_0$epithelial_score, results_0_0_0$epithelial_score),
  group = factor(rep(c("0_1_0", "0_0_0"), 
                     c(nrow(results_0_1_0), nrow(results_0_0_0))))
)
# 
ggplot(combined_data, aes(x = epithelial_score, fill = group)) +
  geom_density(alpha = 0.5, adjust = 0.75) +  # adjust < 1 = less smoothing
  scale_fill_manual(values = c("0_1_0" = "blue", "0_0_0" = "red")) +
  theme_minimal()

# Overlay with transparency
ggplot(combined_data, aes(x = epithelial_score, fill = group)) +
  geom_histogram(bins = 60, alpha = 0.5, position = "identity") +
  scale_fill_manual(values = c("0_1_0" = "blue", "0_0_0" = "red")) +
  theme_minimal() +
  labs(title = "Overlaid Histograms", x = "Epithelial Score", y = "Count")

mean(results_0_1_0$epithelial_score)
mean(results_0_0_0$epithelial_score)
cohens_d(results_0_1_0$epithelial_score, results_0_0_0$epithelial_score)

library(philentropy)

result_fd = calculate_js_divergence(results_0_1_0$epithelial_score, 
                                     results_0_0_0$epithelial_score, 
                                     method =  "fd")

result_fd = calculate_js_divergence(results_0_1_0$epithelial_score, 
                                    results_0_0_0$epithelial_score)
result_fd[1]
result_fd[2]

mean(results_0_0_0$epithelial_score)-mean(results_0_1_0$epithelial_score)
