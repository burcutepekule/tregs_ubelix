rm(list=ls())
library(dplyr)
library(tidyr)
library(zoo)
library(ggplot2)
set.seed(42)

cat("Loading C++ accelerated functions...\n")

cpp_on    = T
one_level = T # not exactly the same result 

param_set_pick = 79902 #19605 VS 79902 - both work
t_max          = 5000
num_reps       = 1 # reps per parameter set

if(cpp_on & one_level){
  source("./MISC/FAST_FUNCTIONS_CPP_ONELEVEL.R")
}else if(cpp_on & !one_level){
 source("./MISC/FAST_FUNCTIONS_CPP.R")
}else{
  source("./MISC/FAST_FUNCTIONS.R")
}

source("./MISC/PLOT_FUNCTIONS.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

cat("\n")

# Optionally check C++ status
# check_cpp_status()  # Uncomment to see which functions are accelerated

# ============================================================================
# SETUP OUTPUT DIRECTORY
# ============================================================================

dir_name_data = './mass_sim_results_R_cpp_macspec_vs_vanilla'
dir.create(dir_name_data, showWarnings = FALSE)

cat("Output directory:", dir_name_data, "\n\n")

colnames_insert = c('epithelial_healthy','epithelial_inj_1','epithelial_inj_2','epithelial_inj_3','epithelial_inj_4','epithelial_inj_5',
                    'phagocyte_M0','phagocyte_M1_L_0','phagocyte_M1_L_1','phagocyte_M1_L_2','phagocyte_M1_L_3','phagocyte_M1_L_4','phagocyte_M1_L_5',
                    'phagocyte_M2_L_0','phagocyte_M2_L_1','phagocyte_M2_L_2','phagocyte_M2_L_3','phagocyte_M2_L_4','phagocyte_M2_L_5',
                    'commensal','pathogen','treg_resting','treg_active','C_ROS','C_M0','C_M1','C_M2','P_ROS','P_M0','P_M1','P_M2')
colnames_insert_1level = c('epithelial_healthy','epithelial_inj_1','epithelial_inj_2','epithelial_inj_3','epithelial_inj_4','epithelial_inj_5',
                    'phagocyte_M0','phagocyte_M1','phagocyte_M2',
                    'commensal','pathogen','treg_resting','treg_active','C_ROS','C_M0','C_M1','C_M2','P_ROS','P_M0','P_M1','P_M2')
# ============================================================================
# READ PARAMETERS FROM CSV
# ============================================================================

cat("Reading parameters...\n")
# params_df = read.csv("./lhs_parameters_ubelix_macspec.csv", stringsAsFactors = FALSE)
# params_df = read.csv("./tregs_better_parameters_ubelix.csv", stringsAsFactors = FALSE)
# params_df = params_df %>% dplyr::filter(param_set_id==1)

params_df = read.csv("./lhs_parameters_ubelix_merged.csv", stringsAsFactors = FALSE)
params_df = params_df %>% dplyr::filter(param_set_id==param_set_pick)
cat("Loaded", nrow(params_df), "parameter sets\n\n")

# if(one_level){ # because it doesn't gradually increase
#   params_df$activity_engulf_M1_baseline = 0.9
#   params_df$activity_engulf_M2_baseline = 0.9
# }

# ============================================================================
# FIXED PARAMETERS (not in CSV)
# ============================================================================

plot_on    = 0
if(plot_on==1){
  dir_name_frames = './frames'
  dir.create(dir_name_frames, showWarnings = FALSE)
}
plot_every = Inf
grid_size  = 25
n_phagocytes = round(grid_size * grid_size * 0.05)
n_tregs = round(grid_size * grid_size * 0.05)
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

cat("Simulation parameters:\n")
cat("  t_max:", t_max, "\n")
cat("  num_reps:", num_reps, "\n")
cat("  grid_size:", grid_size, "x", grid_size, "\n")
cat("  n_phagocytes:", n_phagocytes, "\n")
cat("  n_tregs:", n_tregs, "\n\n")
# ============================================================================
# SCENARIO DEFINITIONS
# ============================================================================

scenarios_df = expand.grid(
  sterile         = c(1),
  allow_tregs     = c(0, 1),
  randomize_tregs = c(0)
)
# DOESN'T MAKE SENSE TO RUN THIS
scenarios_df = scenarios_df %>% dplyr::filter(!(allow_tregs == 0 & randomize_tregs==1))

# ============================================================================
# MAIN SIMULATION LOOP
# ============================================================================
scenario_elapsed_sum = 0
longitudinal_df_keep_all = c()
for(param_set_id_use in 1){
  param_set_use = params_df[param_set_id_use, ]

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
    # RUN SIMULATION WITH C++ ACCELERATION
    # ========================================================================
    
    if(cpp_on & one_level){
      source("./MISC/RUN_REPS_CPP_ONELEVEL.R")
      # colnames(longitudinal_df_keep)[c(9:39)] = colnames_insert_1level 
    }else if(cpp_on & !one_level){
      source("./MISC/RUN_REPS_CPP.R")
      # colnames(longitudinal_df_keep)[c(9:39)] = colnames_insert 
    }else{
      source("./MISC/RUN_REPS.R")
      # colnames(longitudinal_df_keep)[c(9:39)] = colnames_insert 
    }
    
    scenario_end_time = Sys.time()
    scenario_elapsed = as.numeric(difftime(scenario_end_time, scenario_start_time, units = "secs"))


    # saveRDS(longitudinal_df_keep, paste0(dir_name_data,'/longitudinal_df_param_set_id_',param_set_id_use,
    #                                      '_sterile_',sterile,
    #                                      '_tregs_',allow_tregs,
    #                                      '_trnd_',randomize_tregs,'.rds'))

    longitudinal_df_keep_all = rbind(longitudinal_df_keep_all, longitudinal_df_keep)
    cat(sprintf(' - %.1f seconds ✓\n', scenario_elapsed))
    scenario_elapsed_sum = scenario_elapsed_sum + scenario_elapsed
  }
}

variables = c("epithelial_score")

results   = longitudinal_df_keep_all
data_long = results %>%
  dplyr::select(t, sterile, tregs_on, randomize_tregs, rep_id, all_of(variables)) %>%
  pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value")

p = ggplot(data_long, aes(x = t, y = value, color = variable, group = rep_id)) +
  geom_line(alpha = 0.25, linewidth = 1) +
  facet_grid(randomize_tregs ~ sterile + tregs_on , labeller = label_both) +
  scale_color_manual(values = agent_colors) +
  theme_minimal() +
  labs(title = "Epithelial Cell Dynamics", x = "Time", y = "Count", color = "Agent")
print(p)

cat(sprintf('Total: %.1f seconds ✓\n', scenario_elapsed_sum))

