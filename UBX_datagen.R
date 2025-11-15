rm(list=ls())
library(dplyr)
library(tidyr)
# ============================================================================
# LOAD FUNCTIONS
# ============================================================================

source("/storage/homefs/bt25p365/tregs/MISC/FAST_FUNCTIONS.R")
source("/storage/homefs/bt25p365/tregs/MISC/PLOT_FUNCTIONS.R")

split_equal = function(x, n_chunks) {
  split(x, cut(seq_along(x), breaks = n_chunks, labels = FALSE))
}
args   = commandArgs(trailingOnly = TRUE)
n1     = as.integer(args[1])
n2     = as.integer(args[2])

chunks    = split_equal(0:99999, n1)
loop_over = chunks[[n2]]

dir_name_data = '/storage/homefs/bt25p365/tregs/mass_sim_results_R'
dir.create(dir_name_data, showWarnings = FALSE)

colnames_insert = c('epithelial_healthy','epithelial_inj_1','epithelial_inj_2','epithelial_inj_3','epithelial_inj_4','epithelial_inj_5',
                    'phagocyte_M0','phagocyte_M1_L_0','phagocyte_M1_L_1','phagocyte_M1_L_2','phagocyte_M1_L_3','phagocyte_M1_L_4','phagocyte_M1_L_5',
                    'phagocyte_M2_L_0','phagocyte_M2_L_1','phagocyte_M2_L_2','phagocyte_M2_L_3','phagocyte_M2_L_4','phagocyte_M2_L_5',
                    'commensal','pathogen','treg_resting','treg_active','C_ROS','C_M0','C_M1','C_M2','P_ROS','P_M0','P_M1','P_M2')

# ============================================================================
# READ PARAMETERS FROM CSV
# ============================================================================
params_df = read.csv("/storage/homefs/bt25p365/tregs/lhs_parameters_ubelix.csv", stringsAsFactors = FALSE)
params_df = params_df %>% dplyr::filter(param_set_id %in% loop_over)
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

inds_read_filename = paste0('/storage/homefs/bt25p365/tregs/read_ids.rds')
if(file.exists(inds_read_filename)){
  inds_read = readRDS(inds_read_filename)
}else{
  inds_read = c()
}
loop_over = setdiff(loop_over, inds_read)

scenarios_df = expand.grid(
  sterile         = c(0, 1),
  allow_tregs     = c(0, 1),
  randomize_tregs = c(0, 1)
)
# DOESN'T MAKE SENSE TO RUN THIS
scenarios_df = scenarios_df %>% dplyr::filter(!(allow_tregs == 0 & randomize_tregs==1))

for(param_set_id_use in loop_over){
  param_set_use = params_df %>% dplyr::filter(param_set_id==param_set_id_use)
  
  source("/storage/homefs/bt25p365/tregs/MISC/ASSIGN_PARAMETERS.R")
  
  for (scenario_ind in 1:nrow(scenarios_df)){
    sterile         = scenarios_df[scenario_ind,]$sterile
    allow_tregs     = scenarios_df[scenario_ind,]$allow_tregs
    randomize_tregs = scenarios_df[scenario_ind,]$randomize_tregs
    
    print(paste0('Processing param set ',param_set_id_use,' 😱 for scenario ', scenario_ind))
    
    longitudinal_df_keep = c()
    
    source("/storage/homefs/bt25p365/tregs/MISC/RUN_REPS.R")
    
    colnames(longitudinal_df_keep)[c(7:37)] = colnames_insert
    saveRDS(longitudinal_df_keep, paste0(dir_name_data,'/longitudinal_df_param_set_id_',param_set_id_use,
                                         '_sterile_',sterile,
                                         '_tregs_',allow_tregs,
                                         '_trnd_',randomize_tregs,'.rds'))
    
    print(paste0('Data for param set id ',param_set_id_use,' and scenario ', scenario_ind ,' saved 🥳.'))
  }
}




