rm(list=ls())
library(dplyr)
library(tidyr)
library(zoo)
library(cowplot)
library(av)
library(ggplot2)

setwd('~/Dropbox/tregs_ubelix')
source("./MISC/FAST_FUNCTIONS.R")
source("./MISC/PLOT_FUNCTIONS.R")
source("./MISC/DATA_READ_FUNCTIONS.R")

# params_df    = read.csv("./tregs_better_parameters_ubelix.csv", stringsAsFactors = FALSE)
# loop_over    = c(3000)
params_df    = read.csv("./lhs_parameters_ubelix.csv", stringsAsFactors = FALSE)
loop_over    = c(5)
params_df    = params_df %>% dplyr::filter(param_set_id %in% loop_over)

# ============================================================================
# SETUP OUTPUT DIRECTORY
# ============================================================================

dir_name_data = '/Users/burcutepekule/Desktop/gif_out'
dir.create(dir_name_data, showWarnings = FALSE)

cat("Output directory:", dir_name_data, "\n\n")

colnames_insert = c('epithelial_healthy','epithelial_inj_1','epithelial_inj_2','epithelial_inj_3','epithelial_inj_4','epithelial_inj_5',
                    'phagocyte_M0','phagocyte_M1_L_0','phagocyte_M1_L_1','phagocyte_M1_L_2','phagocyte_M1_L_3','phagocyte_M1_L_4','phagocyte_M1_L_5',
                    'phagocyte_M2_L_0','phagocyte_M2_L_1','phagocyte_M2_L_2','phagocyte_M2_L_3','phagocyte_M2_L_4','phagocyte_M2_L_5',
                    'commensal','pathogen','treg_resting','treg_active','C_ROS','C_M0','C_M1','C_M2','P_ROS','P_M0','P_M1','P_M2')

# ============================================================================
# FIXED PARAMETERS (not in CSV)
# ============================================================================
num_reps   = 1
t_max      = 1000
plot_on    = 1
if(plot_on==1){
  dir_name_frames = paste0(dir_name_data,'/frames')
  dir.create(dir_name_frames, showWarnings = FALSE)
}
plot_every = 100
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

# === All you need for smart macrophages
scenarios_df = expand.grid(
  sterile         = c(0), 
  allow_tregs     = c(0, 1),
  randomize_tregs = c(0, 1)
)
scenarios_df = scenarios_df %>% dplyr::filter(!(allow_tregs==0 & randomize_tregs==1))

cat("Running", nrow(scenarios_df), "scenarios per parameter set\n")
cat("Total simulations:", length(loop_over) * nrow(scenarios_df) * num_reps, "\n\n")

# ============================================================================
# MAIN SIMULATION LOOP
# ============================================================================

for(param_set_id_use in loop_over){
  param_set_use = params_df %>% dplyr::filter(param_set_id==param_set_id_use)

  for (scenario_ind in 2:nrow(scenarios_df)){
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

    
    # Define file list and output
    pattern   = paste0("^frame_param_",param_set_id_use,"_rep_0_STERILE_", sterile, "_TREGS_", allow_tregs, "_trnd_", randomize_tregs, "_\\d+\\.png$")
    png_files = list.files(paste0(dir_name_frames,"/"), full.names = TRUE, pattern = pattern)
    # Sort files numerically by the time value (last number before .png)
    png_files = png_files[order(as.numeric(gsub(".*_(\\d+)\\.png$", "\\1", png_files)))]
    video_out = paste0(dir_name_frames,"/simulation_param_",param_set_id_use,"_sterile", sterile, "_tregs_", allow_tregs, "_trnd_", randomize_tregs, ".mp4")
    
    # Create video
    av_encode_video(
      input = png_files,
      output = video_out,
      framerate = 3, # slower playback (e.g. 2 FPS)
      vfilter = "scale=1000:-2",  # Resize if needed
      codec = "libx264"      # H.264 codec is widely supported
    )
    
    cat(sprintf(' - %.1f seconds ✓\n', scenario_elapsed))
  }
}

