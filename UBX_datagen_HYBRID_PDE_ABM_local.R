rm(list=ls())
library(dplyr)
library(tidyr)
library(zoo)
set.seed(42)

cat("Loading C++ accelerated functions...\n")
source("./tregs/MISC/FAST_FUNCTIONS_CPP_ONELEVEL.R")
source("./tregs/MISC/PLOT_FUNCTIONS.R")
source("./tregs/MISC/DATA_READ_FUNCTIONS.R")

cat("\n")
cat("=" , rep("=", 70), sep = "")
cat("\n")
cat("HYBRID PDE-ABM MODEL\n")
cat("Microbes: PDE (continuous concentration fields)\n")
cat("Agents: ABM (discrete epithelial cells, macrophages, Tregs)\n")
cat("Memory: Exponential decay (biologically realistic)\n")
cat("=" , rep("=", 70), sep = "")
cat("\n\n")

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

dir_name_data = './tregs/mass_sim_results_HYBRID_PDE_ABM'
dir.create(dir_name_data, showWarnings = FALSE)

cat("Output directory:", dir_name_data, "\n\n")

colnames_insert_1level = c('epithelial_healthy','epithelial_inj_1','epithelial_inj_2',
                           'epithelial_inj_3','epithelial_inj_4','epithelial_inj_5',
                           'phagocyte_M0','phagocyte_M1','phagocyte_M2',
                           'commensal','pathogen','treg_resting','treg_active',
                           'C_ROS','C_M0','C_M1','C_M2','P_ROS','P_M0','P_M1','P_M2')

# ============================================================================
# READ PARAMETERS FROM CSV
# ============================================================================

cat("Reading parameters...\n")
params_df = read.csv("./tregs/lhs_parameters_ubelix_macspec.csv", stringsAsFactors = FALSE)
params_df = params_df %>% dplyr::filter(param_set_id %in% loop_over)
cat("Loaded", nrow(params_df), "parameter sets\n\n")

# ============================================================================
# FIXED PARAMETERS (not in CSV)
# ============================================================================
t_max      = 5000
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

# NEW HYBRID MODEL PARAMETERS
memory_half_life = 2  # timesteps for engulfment memory decay
D_microbe        = 0.0625 # diffusion rate for pathogens/commensals -> 0.5/8
max_microbe_concentration = 100  # maximum microbe concentration per cell

cat("Simulation parameters:\n")
cat("  t_max:", t_max, "\n")
cat("  grid_size:", grid_size, "x", grid_size, "\n")
cat("  n_phagocytes:", n_phagocytes, "\n")
cat("  n_tregs:", n_tregs, "\n")
cat("  memory_half_life:", memory_half_life, "timesteps\n")
cat("  D_microbe:", D_microbe, "\n\n")

# ============================================================================
# SCENARIO DEFINITIONS
# ============================================================================

scenarios_df = expand.grid(
  control         = c(0),
  sterile         = c(0, 1),
  allow_tregs     = c(0, 1),
  randomize_tregs = c(0, 1),
  macspec_on      = c(0, 1)
)
scenarios_df = scenarios_df %>% dplyr::filter(!(allow_tregs == 0 & randomize_tregs==1))
scenarios_df = scenarios_df %>% dplyr::filter(!(macspec_on ==1 & allow_tregs == 1 & randomize_tregs==1))
scenarios_df = scenarios_df %>% dplyr::filter(!(macspec_on ==1 & allow_tregs == 1 & randomize_tregs==0))

scenarios_df_ctrl = expand.grid(
  control         = c(1),
  sterile         = c(0, 1),
  allow_tregs     = c(0),
  randomize_tregs = c(0),
  macspec_on      = c(0)
)
scenarios_df=rbind(scenarios_df_ctrl, scenarios_df)

cat("Total scenarios:", nrow(scenarios_df), "\n\n")

# ============================================================================
# MAIN SIMULATION LOOP
# ============================================================================

for(param_set_id_use in loop_over){
  param_set_use = params_df %>% dplyr::filter(param_set_id==param_set_id_use)

  for (scenario_ind in 1:nrow(scenarios_df)){
    sterile         = scenarios_df[scenario_ind,]$sterile
    allow_tregs     = scenarios_df[scenario_ind,]$allow_tregs
    randomize_tregs = scenarios_df[scenario_ind,]$randomize_tregs
    macspec_on      = scenarios_df[scenario_ind,]$macspec_on
    control         = scenarios_df[scenario_ind,]$control

    # if(control==1){
    #   num_reps = 1
    # }else{
    #   num_reps = 100
    # }
    num_reps = 100
    
    source("./MISC/ASSIGN_PARAMETERS_MACSPEC.R")
    cat("  num_reps:", num_reps, "\n")

    cat(paste0('[', Sys.time(), '] Processing param set ', param_set_id_use,
               ' - scenario ', scenario_ind, '/', nrow(scenarios_df)))

    scenario_start_time  = Sys.time()

    longitudinal_df_keep = c()

    # ========================================================================
    # RUN HYBRID PDE-ABM SIMULATION
    # ========================================================================

    source("./MISC/RUN_REPS_HYBRID_PDE_ABM.R")

    scenario_end_time = Sys.time()
    scenario_elapsed  = as.numeric(difftime(scenario_end_time, scenario_start_time, units = "secs"))

    saveRDS(longitudinal_df_keep, paste0(dir_name_data,'/longitudinal_df_param_set_id_',param_set_id_use,
                                         '_control_',control,
                                         '_sterile_',sterile,
                                         '_macspec_',macspec_on,
                                         '_tregs_',allow_tregs,
                                         '_trnd_',randomize_tregs,'.rds'))

    cat(sprintf(' - %.1f seconds ✓\n', scenario_elapsed))
  }
}

# ============================================================================
# SUMMARY
# ============================================================================
cat("\n")
cat("=" , rep("=", 70), sep = "")
cat("\n")
cat("SIMULATION COMPLETE\n")
cat("=" , rep("=", 70), sep = "")
cat("\n")
cat("Data saved to:", dir_name_data, "\n")
cat("\n")
