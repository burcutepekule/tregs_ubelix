rm(list=ls())
library(dplyr)
library(tidyr)
library(zoo)

source("./MISC/FAST_FUNCTIONS.R")
source("./MISC/PLOT_FUNCTIONS.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

# ============================================================================
# loop_over = c(81709) # param IDs of choice
loop_over = c(81709) # param IDs of choice

# ============================================================================

dir_name_data = './mass_sim_results_R_cpp_macspec_local'
dir.create(dir_name_data, showWarnings = FALSE)

cat("Output directory:", dir_name_data, "\n\n")

colnames_insert = c('epithelial_healthy','epithelial_inj_1','epithelial_inj_2','epithelial_inj_3','epithelial_inj_4','epithelial_inj_5',
                    'phagocyte_M0','phagocyte_M1_L_0','phagocyte_M1_L_1','phagocyte_M1_L_2','phagocyte_M1_L_3','phagocyte_M1_L_4','phagocyte_M1_L_5',
                    'phagocyte_M2_L_0','phagocyte_M2_L_1','phagocyte_M2_L_2','phagocyte_M2_L_3','phagocyte_M2_L_4','phagocyte_M2_L_5',
                    'commensal','pathogen','treg_resting','treg_active','C_ROS','C_M0','C_M1','C_M2','P_ROS','P_M0','P_M1','P_M2')

params_df = read.csv("./lhs_parameters_ubelix.csv", stringsAsFactors = FALSE)
params_df = params_df %>% dplyr::filter(param_set_id %in% loop_over)
cat("New parameters: mac_discrimination_efficiency, mac_rat_com_pat_threshold\n\n")

# ============================================================================
# FIXED PARAMETERS (not in CSV)
# ============================================================================

t_max      = 1000
num_reps   = 1 # reps per parameter set
plot_on    = 0
if(plot_on==1){
  dir_name_frames = './frames'
  dir.create(dir_name_frames, showWarnings = FALSE)
}
plot_every = Inf
grid_size  = 25
n_phagocytes = round(grid_size * grid_size * 0.35)
n_tregs = round(grid_size * grid_size * 0.35)
n_commensals_lp = 20

injury_percentage = 60
max_level_injury  = 5

max_cell_value_ROS   = 1
max_cell_value_DAMPs = 1
max_cell_value_SAMPs = 1

lim_ROS  = max_cell_value_ROS
lim_DAMP = max_cell_value_DAMPs
lim_SAMP = max_cell_value_SAMPs

act_radius_ROS   = 1
act_radius_treg  = 1
act_radius_DAMPs = 1
act_radius_SAMPs = 1

# Logistic function parameters (for epithelial injury calculation)
k_in  = 0.044
x0_in = 50

# ============================================================================
# SCENARIO DEFINITIONS
# ============================================================================

scenarios_df = expand.grid(
  sterile         = c(0, 1),
  allow_tregs     = c(0, 1),
  randomize_tregs = c(0, 1)
)
# DOESN'T MAKE SENSE TO RUN THIS
scenarios_df = scenarios_df %>% dplyr::filter(!(allow_tregs == 0 & randomize_tregs==1))

# ============================================================================
# MAIN SIMULATION LOOP
# ============================================================================
longitudinal_df_keep_all = c()
for(param_set_id_use in loop_over){
  param_set_use = params_df %>% dplyr::filter(param_set_id==param_set_id_use)
  
  for (scenario_ind in 1:nrow(scenarios_df)){
    sterile         = scenarios_df[scenario_ind,]$sterile
    allow_tregs     = scenarios_df[scenario_ind,]$allow_tregs
    randomize_tregs = scenarios_df[scenario_ind,]$randomize_tregs
    
    source("./MISC/ASSIGN_PARAMETERS.R")
    
    cat(paste0('[', Sys.time(), '] Processing param set ', param_set_id_use,
               ' - scenario ', scenario_ind, '/', nrow(scenarios_df)))
    
    # Track timing for this scenario
    scenario_start_time = Sys.time()
    
    longitudinal_df_keep = c()
    
    # ========================================================================
    # RUN SIMULATION WITH C++ ACCELERATION AND MACROPHAGE SPECIFICITY
    # ========================================================================
    source("./MISC/RUN_REPS_CPP.R")
    
    scenario_end_time = Sys.time()
    scenario_elapsed = as.numeric(difftime(scenario_end_time, scenario_start_time, units = "secs"))
    
    colnames(longitudinal_df_keep)[c(9:39)] = colnames_insert # this sould be 9:39 - fix it when reading
    # saveRDS(longitudinal_df_keep, paste0(dir_name_data,'/longitudinal_df_param_set_id_',param_set_id_use,
    #                                      '_sterile_',sterile,
    #                                      '_tregs_',allow_tregs,
    #                                      '_trnd_',randomize_tregs,'.rds'))
    longitudinal_df_keep_all = rbind(longitudinal_df_keep_all, longitudinal_df_keep)
    
    cat(sprintf(' - %.1f seconds ✓\n', scenario_elapsed))
  }
}
results_vanilla = longitudinal_df_keep_all

variables = c("epithelial_score")
data_long_vanilla = results_vanilla %>%
  dplyr::select(t, sterile, tregs_on, randomize_tregs, rep_id, all_of(variables)) %>%
  pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value")

p_vanilla = ggplot(data_long_vanilla, aes(x = t, y = value, color = variable, group = rep_id)) +
  geom_line(alpha = 0.5, linewidth = 1) +
  facet_grid(randomize_tregs ~ sterile + tregs_on , labeller = label_both) +
  scale_color_manual(values = agent_colors) +
  theme_minimal() +
  labs(title = "Epithelial Cell Dynamics", x = "Time", y = "Count", color = "Agent")

# ========================================================================
scenarios_df = scenarios_df %>% dplyr::filter(allow_tregs == 0) #the point is NOT to have Tregs!

longitudinal_df_keep_all = c()
for(param_set_id_use in loop_over){
  param_set_use = params_df %>% dplyr::filter(param_set_id==param_set_id_use)

  # ================= These parameters are not in the original dataset =========
  param_set_use$mac_discrimination_efficiency = 1
  param_set_use$mac_rat_com_pat_threshold     = 0.75
  # ================= These parameters are not in the original dataset =========
  
  for (scenario_ind in 1:nrow(scenarios_df)){
    sterile         = scenarios_df[scenario_ind,]$sterile
    allow_tregs     = scenarios_df[scenario_ind,]$allow_tregs
    randomize_tregs = scenarios_df[scenario_ind,]$randomize_tregs

    source("./MISC/ASSIGN_PARAMETERS_MACSPEC.R")

    cat(paste0('[', Sys.time(), '] Processing param set ', param_set_id_use,
               ' - scenario ', scenario_ind, '/', nrow(scenarios_df)))

    # Track timing for this scenario
    scenario_start_time = Sys.time()

    longitudinal_df_keep = c()

    # ========================================================================
    # RUN SIMULATION WITH C++ ACCELERATION AND MACROPHAGE SPECIFICITY
    # ========================================================================
    source("./MISC/RUN_REPS_CPP_MACSPEC.R")

    scenario_end_time = Sys.time()
    scenario_elapsed = as.numeric(difftime(scenario_end_time, scenario_start_time, units = "secs"))

    colnames(longitudinal_df_keep)[c(9:39)] = colnames_insert 
    # saveRDS(longitudinal_df_keep, paste0(dir_name_data,'/longitudinal_df_param_set_id_',param_set_id_use,
    #                                      '_sterile_',sterile,
    #                                      '_tregs_',allow_tregs,
    #                                      '_trnd_',randomize_tregs,'.rds'))
    longitudinal_df_keep_all = rbind(longitudinal_df_keep_all, longitudinal_df_keep)
    
    cat(sprintf(' - %.1f seconds ✓\n', scenario_elapsed))
  }
}
results_macspec = longitudinal_df_keep_all

variables = c("epithelial_score")
data_long_macspec = results_macspec %>%
  dplyr::select(t, sterile, tregs_on, randomize_tregs, rep_id, all_of(variables)) %>%
  pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value")

p_macspec = ggplot(data_long_macspec, aes(x = t, y = value, color = variable, group = rep_id)) +
  geom_line(alpha = 0.5, linewidth = 1) +
  facet_grid(randomize_tregs ~ sterile + tregs_on , labeller = label_both) +
  scale_color_manual(values = agent_colors) +
  theme_minimal() +
  labs(title = "Epithelial Cell Dynamics", x = "Time", y = "Count", color = "Agent")

#========== PLOTTING ===========================================================
data_long_vanilla$type = 'Vanilla'
data_long_macspec$type = 'MACSPEC'
data_long = rbind(data_long_vanilla, data_long_macspec)

p_merged = ggplot(data_long, aes(x = t, y = value, color = type, group = interaction(type, rep_id))) +
  geom_line(alpha = 0.5, linewidth = 1) +
  scale_color_manual(values = c("Vanilla" = "#B6F500", "MACSPEC" = "#4300FF")) +
  facet_grid(randomize_tregs ~ sterile + tregs_on , labeller = label_both) +
  theme_minimal() +
  labs(title = "Epithelial Cell Dynamics", x = "Time", y = "Count", color = "Agent")
print(p_merged)
