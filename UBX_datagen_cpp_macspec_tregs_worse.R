rm(list=ls())
library(dplyr)
library(tidyr)
library(zoo)
setwd('~/Dropbox/tregs_ubelix')
source("./MISC/FAST_FUNCTIONS.R")
source("./MISC/PLOT_FUNCTIONS.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

split_equal = function(x, n_chunks) {
  split(x, cut(seq_along(x), breaks = n_chunks, labels = FALSE))
}

# ============================================================================
# COMMAND LINE ARGUMENTS
# ============================================================================

args   = commandArgs(trailingOnly = TRUE)
n1     = as.integer(args[1])
n2     = as.integer(args[2])

# ============================================================================
# READ PARAMETERS FROM CSV
# ============================================================================

cat("Reading parameters with macrophage specificity for tregs better...\n")
params_df = read.csv("./lhs_parameters_macspec_expanded_tregs_worse.csv", stringsAsFactors = FALSE)

# ============================================================================
# SPLIT
# ============================================================================

chunks    = split_equal(params_df$param_set_id_new, n1)
loop_over = chunks[[n2]]
params_df = params_df %>% dplyr::filter(param_set_id_new %in% loop_over)

cat("Processing chunk", n2, "of", n1, "\n")
cat("Parameter sets:", loop_over[1], "-", loop_over[length(loop_over)], "\n\n")

# ============================================================================
# SETUP OUTPUT DIRECTORY
# ============================================================================

dir_name_data = '/Users/burcutepekule/Desktop/mass_sim_results_R_cpp_macspec_tregs_worse'
dir.create(dir_name_data, showWarnings = FALSE)

cat("Output directory:", dir_name_data, "\n\n")

colnames_insert = c('epithelial_healthy','epithelial_inj_1','epithelial_inj_2','epithelial_inj_3','epithelial_inj_4','epithelial_inj_5',
                    'phagocyte_M0','phagocyte_M1_L_0','phagocyte_M1_L_1','phagocyte_M1_L_2','phagocyte_M1_L_3','phagocyte_M1_L_4','phagocyte_M1_L_5',
                    'phagocyte_M2_L_0','phagocyte_M2_L_1','phagocyte_M2_L_2','phagocyte_M2_L_3','phagocyte_M2_L_4','phagocyte_M2_L_5',
                    'commensal','pathogen','treg_resting','treg_active','C_ROS','C_M0','C_M1','C_M2','P_ROS','P_M0','P_M1','P_M2')



# ============================================================================
# FIXED PARAMETERS (not in CSV)
# ============================================================================

t_max      = 5000
num_reps   = 10 # reps per parameter set
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

cat("Simulation parameters:\n")
cat("  t_max:", t_max, "\n")
cat("  num_reps:", num_reps, "\n")
cat("  grid_size:", grid_size, "x", grid_size, "\n")
cat("  n_phagocytes:", n_phagocytes, "\n")
cat("  n_tregs:", n_tregs, "\n\n")

# ============================================================================
# SCENARIO DEFINITIONS
# ============================================================================

# === All you need for smart macrophages
scenarios_df = expand.grid(
  sterile = c(0, 1),
  allow_tregs     = c(0),
  randomize_tregs = c(0)
)

cat("Running", nrow(scenarios_df), "scenarios per parameter set\n")
cat("Total simulations:", length(loop_over) * nrow(scenarios_df) * num_reps, "\n\n")

# ============================================================================
# MAIN SIMULATION LOOP
# ============================================================================

for(param_set_id_use in loop_over){
  param_set_use = params_df %>% dplyr::filter(param_set_id_new==param_set_id_use)

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

    colnames(longitudinal_df_keep)[c(9:39)] = colnames_insert # this sould be 9:39 - fix it when reading
    saveRDS(longitudinal_df_keep, paste0(dir_name_data,'/longitudinal_df_param_set_id_',param_set_id_use,
                                         '_sterile_',sterile,
                                         '_tregs_',allow_tregs,
                                         '_trnd_',randomize_tregs,'.rds'))

    cat(sprintf(' - %.1f seconds ✓\n', scenario_elapsed))
  }
}

