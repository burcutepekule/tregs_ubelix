rm(list=ls())
library(dplyr)
library(tidyr)

# ============================================================================
# C++ ACCELERATED DATA GENERATION SCRIPT
# ============================================================================
# This script uses C++ implementations for 20-100x speedup
# To use: Rscript UBX_datagen_cpp.R <n1> <n2>
# ============================================================================

cat("\n")
cat("=" = rep("=", 70), sep = "")
cat("\n🚀 C++ ACCELERATED TREGS SIMULATION\n")
cat("=" = rep("=", 70), sep = "")
cat("\n\n")

# ============================================================================
# LOAD C++ ACCELERATED FUNCTIONS
# ============================================================================

cat("Loading C++ accelerated functions...\n")
source("MISC/FAST_FUNCTIONS_CPP.R")
source("MISC/PLOT_FUNCTIONS.R")

cat("\n")

# Optionally check C++ status
# check_cpp_status()  # Uncomment to see which functions are accelerated

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

chunks    = split_equal(0:99999, n1)
loop_over = chunks[[n2]]

cat("Processing chunk", n2, "of", n1, "\n")
cat("Parameter sets:", min(loop_over), "-", max(loop_over), "\n\n")

# ============================================================================
# SETUP OUTPUT DIRECTORY
# ============================================================================

dir_name_data = '/storage/homefs/bt25p365/tregs/mass_sim_results_R_cpp'
dir.create(dir_name_data, showWarnings = FALSE)

cat("Output directory:", dir_name_data, "\n\n")

colnames_insert = c('epithelial_healthy','epithelial_inj_1','epithelial_inj_2','epithelial_inj_3','epithelial_inj_4','epithelial_inj_5',
                    'phagocyte_M0','phagocyte_M1_L_0','phagocyte_M1_L_1','phagocyte_M1_L_2','phagocyte_M1_L_3','phagocyte_M1_L_4','phagocyte_M1_L_5',
                    'phagocyte_M2_L_0','phagocyte_M2_L_1','phagocyte_M2_L_2','phagocyte_M2_L_3','phagocyte_M2_L_4','phagocyte_M2_L_5',
                    'commensal','pathogen','treg_resting','treg_active','C_ROS','C_M0','C_M1','C_M2','P_ROS','P_M0','P_M1','P_M2')

# ============================================================================
# READ PARAMETERS FROM CSV
# ============================================================================

cat("Reading parameters...\n")
params_df = read.csv("lhs_parameters_ubelix.csv", stringsAsFactors = FALSE)
params_df = params_df %>% dplyr::filter(param_set_id %in% loop_over)
cat("Loaded", nrow(params_df), "parameter sets\n\n")

# ============================================================================
# FIXED PARAMETERS (not in CSV)
# ============================================================================

t_max      = 5000
num_reps   = 100 # reps per parameter set
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
# READ TRACKING FILE
# ============================================================================

inds_read_filename = paste0('/storage/homefs/bt25p365/tregs/read_ids.rds')
if(file.exists(inds_read_filename)){
  inds_read = readRDS(inds_read_filename)
}else{
  inds_read = c()
}
loop_over = setdiff(loop_over, inds_read)

cat("After filtering already processed:", length(loop_over), "parameter sets remaining\n\n")

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

cat("Running", nrow(scenarios_df), "scenarios per parameter set\n")
cat("Total simulations:", length(loop_over) * nrow(scenarios_df) * num_reps, "\n\n")

cat("=" = rep("=", 70), sep = "")
cat("\n")
cat("Starting simulations with C++ acceleration...\n")
cat("=" = rep("=", 70), sep = "")
cat("\n\n")

# ============================================================================
# MAIN SIMULATION LOOP
# ============================================================================

start_time_total = Sys.time()

for(param_set_id_use in loop_over){
  param_set_use = params_df %>% dplyr::filter(param_set_id==param_set_id_use)

  source("MISC/ASSIGN_PARAMETERS.R")

  for (scenario_ind in 1:nrow(scenarios_df)){
    sterile         = scenarios_df[scenario_ind,]$sterile
    allow_tregs     = scenarios_df[scenario_ind,]$allow_tregs
    randomize_tregs = scenarios_df[scenario_ind,]$randomize_tregs

    cat(paste0('[', Sys.time(), '] Processing param set ', param_set_id_use,
               ' - scenario ', scenario_ind, '/', nrow(scenarios_df)))

    # Track timing for this scenario
    scenario_start_time = Sys.time()

    longitudinal_df_keep = c()

    # ========================================================================
    # RUN SIMULATION WITH C++ ACCELERATION
    # ========================================================================
    source("MISC/RUN_REPS_CPP.R")

    scenario_end_time = Sys.time()
    scenario_elapsed = as.numeric(difftime(scenario_end_time, scenario_start_time, units = "secs"))

    colnames(longitudinal_df_keep)[c(7:37)] = colnames_insert
    saveRDS(longitudinal_df_keep, paste0(dir_name_data,'/longitudinal_df_param_set_id_',param_set_id_use,
                                         '_sterile_',sterile,
                                         '_tregs_',allow_tregs,
                                         '_trnd_',randomize_tregs,'.rds'))

    cat(sprintf(' - %.1f seconds ✓\n', scenario_elapsed))
  }
}

# ============================================================================
# SUMMARY
# ============================================================================

end_time_total = Sys.time()
total_elapsed = as.numeric(difftime(end_time_total, start_time_total, units = "secs"))

cat("\n")
cat("=" = rep("=", 70), sep = "")
cat("\n")
cat("SIMULATION COMPLETE\n")
cat("=" = rep("=", 70), sep = "")
cat("\n")
cat("Total time:", sprintf("%.1f seconds (%.1f minutes)\n", total_elapsed, total_elapsed/60))
cat("Data saved to:", dir_name_data, "\n")
cat("\n")

cat("🚀 C++ acceleration enabled!\n")
cat("   Expected speedup: 20-100x faster than pure R\n")
cat("=" = rep("=", 70), sep = "")
cat("\n\n")
