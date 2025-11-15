rm(list=ls())
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

source("/Users/burcutepekule/Dropbox/tregs_ubelix/MISC/PLOT_FUNCTIONS.R")

path         = "/Users/burcutepekule/Desktop/mass_sim_results_R_cpp/"
param_id     = 0 
rep_ind      = 9

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

# Merge
results = results %>% dplyr::mutate(epithelial_score = 6*epithelial_healthy+ # higher the score, healthier the epithelium!
                                      5*epithelial_inj_1+
                                      4*epithelial_inj_2+
                                      3*epithelial_inj_3+
                                      2*epithelial_inj_4+
                                      1*epithelial_inj_5)

results = results %>% dplyr::filter(rep_id==rep_ind)
variables = c("epithelial_healthy", paste0("epithelial_inj_", 1:5))
data_long = results %>%
  dplyr::select(t, sterile, tregs_on, randomize_tregs, all_of(variables)) %>%
  pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value")

p=ggplot(data_long, aes(x = t, y = value, color = variable)) +
  geom_line(alpha = 1, linewidth = 1) +
  facet_grid(randomize_tregs ~ sterile + tregs_on , labeller = label_both) +
  scale_color_manual(values = agent_colors) +
  theme_minimal() +
  labs(title = "Epithelial Cell Dynamics", x = "Time", y = "Count", color = "Agent")
print(p)



