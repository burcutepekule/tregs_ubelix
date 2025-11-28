rm(list=ls())
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

source("/Users/burcutepekule/Dropbox/tregs_ubelix/MISC/PLOT_FUNCTIONS.R")


colnames_insert_fixed = c("t","sterile","tregs_on","randomize_tregs","param_set_id","rep_id","epithelial_score","time_ss",
                          'epithelial_healthy','epithelial_inj_1','epithelial_inj_2','epithelial_inj_3','epithelial_inj_4',
                          'epithelial_inj_5','phagocyte_M0','phagocyte_M1_L_0','phagocyte_M1_L_1','phagocyte_M1_L_2','phagocyte_M1_L_3',
                          'phagocyte_M1_L_4','phagocyte_M1_L_5','phagocyte_M2_L_0','phagocyte_M2_L_1','phagocyte_M2_L_2',
                          'phagocyte_M2_L_3','phagocyte_M2_L_4','phagocyte_M2_L_5','commensal','pathogen','treg_resting',
                          'treg_active','C_ROS','C_M0','C_M1','C_M2','P_ROS','P_M0','P_M1','P_M2')

path         = "/Users/burcutepekule/Desktop/mass_sim_results_R_cpp/"
param_id     = 29510

for(rep_ind_vec in 0:99){
  
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
  
  results_000 = results %>% dplyr::filter(sterile==0 & tregs_on==0 & randomize_tregs==0)
  results_010 = results %>% dplyr::filter(sterile==0 & tregs_on==1 & randomize_tregs==0)
  
  results_keep = results
  
  if(!any(is.na(rep_ind_vec))){
    results = results_keep %>% dplyr::filter(rep_id %in% rep_ind_vec)  
  }
  
  # variables = c("epithelial_healthy", paste0("epithelial_inj_", 1:5))
  variables = c("epithelial_score")
  
  data_long = results %>%
    dplyr::select(t, sterile, tregs_on, randomize_tregs, rep_id, all_of(variables)) %>%
    pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value")
  
  
  browser()
  
  p=ggplot(data_long, aes(x = t, y = value, color = variable, group = rep_id)) +
    geom_line(alpha = 1, linewidth = 1) +
    facet_grid(randomize_tregs ~ sterile + tregs_on , labeller = label_both) +
    scale_color_manual(values = agent_colors) +
    theme_minimal() +
    labs(title = "Epithelial Cell Dynamics", x = "Time", y = "Count", color = "Agent")
  
  print(p)
  print(rep_ind_vec)
  
}
