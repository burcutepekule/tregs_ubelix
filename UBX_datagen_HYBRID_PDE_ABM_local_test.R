rm(list=ls())
library(dplyr)
library(tidyr)
library(zoo)
set.seed(42)

cat("Loading C++ accelerated functions...\n")
source("./MISC/FAST_FUNCTIONS_CPP_ONELEVEL.R")
source("./MISC/PLOT_FUNCTIONS.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

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
# SETUP OUTPUT DIRECTORY
# ============================================================================

dir_name_data = './mass_sim_results_HYBRID_PDE_ABM'
dir.create(dir_name_data, showWarnings = FALSE)

cat("Output directory:", dir_name_data, "\n\n")

colnames_insert_1level = c('epithelial_healthy','epithelial_inj_1','epithelial_inj_2','epithelial_inj_3','epithelial_inj_4','epithelial_inj_5',
                           'phagocyte_M0','phagocyte_M1','phagocyte_M2',
                           'commensal','pathogen','treg_resting','treg_active','C_ROS','C_M0','C_M1','C_M2','P_ROS','P_M0','P_M1','P_M2')

# ============================================================================
# READ PARAMETERS FROM CSV
# ============================================================================

cat("Reading parameters...\n")
loop_over = c(1)
params_df = read.csv("./lhs_parameters_ubelix_macspec.csv", stringsAsFactors = FALSE)
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
memory_half_life = 5  # timesteps for engulfment memory decay
D_microbe = 0.1      # diffusion rate for pathogens/commensals (max 0.12!)
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
scenario_elapsed_sum = 0
longitudinal_df_keep_all = c()

for(param_set_id_use in loop_over){
  param_set_use = params_df %>% dplyr::filter(param_set_id==param_set_id_use)

  for (scenario_ind in 1:nrow(scenarios_df)){
    sterile         = scenarios_df[scenario_ind,]$sterile
    allow_tregs     = scenarios_df[scenario_ind,]$allow_tregs
    randomize_tregs = scenarios_df[scenario_ind,]$randomize_tregs
    macspec_on      = scenarios_df[scenario_ind,]$macspec_on
    control         = scenarios_df[scenario_ind,]$control

    if(control==1){
      num_reps = 1
    }else{
      num_reps = 100
    }

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
    
    longitudinal_df_keep_all = rbind(longitudinal_df_keep_all, longitudinal_df_keep)
    scenario_elapsed_sum = scenario_elapsed_sum + scenario_elapsed
  }
}

variables = c("pathogen")

results   = longitudinal_df_keep_all
data_long = results %>%
  dplyr::select(t, control, sterile, tregs_on, macspec_on, randomize_tregs, rep_id, all_of(variables)) %>%
  pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value")

p_p = ggplot(data_long, aes(x = t, y = value, color = variable, group = rep_id)) +
  geom_line(alpha = 0.25, linewidth = 1) +
  facet_grid(randomize_tregs ~ control + sterile + macspec_on + tregs_on , labeller = label_both) +
  scale_color_manual(values = agent_colors) +
  theme_minimal() +
  labs(title = "Pathogens", x = "Time", y = "Count", color = "Agent")

variables = c("epithelial_score")

results   = longitudinal_df_keep_all
data_long = results %>%
  dplyr::select(t, control, sterile, tregs_on, macspec_on, randomize_tregs, rep_id, all_of(variables)) %>%
  pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value")

p_e = ggplot(data_long, aes(x = t, y = value, color = variable, group = rep_id)) +
  geom_line(alpha = 0.25, linewidth = 1) +
  facet_grid(randomize_tregs ~ control + sterile + macspec_on + tregs_on , labeller = label_both) +
  scale_color_manual(values = agent_colors) +
  theme_minimal() +
  labs(title = "Epithelial Cells", x = "Time", y = "Count", color = "Agent")

cowplot::plot_grid(
  p_p, p_e,      # your two plots
  ncol = 1,      # 1 column
  nrow = 2,      # 2 rows (optional, inferred from ncol + number of plots)
  align = "v",   # vertically align
  rel_heights = c(1, 1)  # adjust if you want one taller than the other
)

cat(sprintf('Total: %.1f seconds ✓\n', scenario_elapsed_sum))


