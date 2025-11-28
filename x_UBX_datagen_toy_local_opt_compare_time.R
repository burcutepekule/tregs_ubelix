rm(list=ls())
library(dplyr)
library(tidyr)

setwd('/Users/burcutepekule/Dropbox/tregs_ubelix')
source("./MISC/FAST_FUNCTIONS.R")


t_max      = 5000      # Increase this for more realistic timing (100-500)
num_reps   = 1       # Number of replicates per version
plot_on    = 0       # Keep plotting off for timing
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
# Logistic function parameters
k_in  = 0.044
x0_in = 50


# Use a single parameter set for controlled comparison

sterile         = 0
allow_tregs     = 1
randomize_tregs = 1
# Load a test parameter set (or define manually)

params_df = read.csv("./lhs_parameters_local.csv", stringsAsFactors = FALSE)
param_set_id_use = 2
param_set_use = params_df %>% dplyr::filter(param_set_id == param_set_id_use)

# Assign parameters

source("MISC/ASSIGN_PARAMETERS.R")

# ============================================================================
# RUN ORIGINAL VERSION WITH TIMING
# ============================================================================

cat("\n")
cat("=" = rep("=", 70), sep="")
cat("\nRUNNING ORIGINAL VERSION (RUN_REPS.R)\n")
cat("=" = rep("=", 70), sep="")
cat("\n")
cat("Parameters: t_max =", t_max, ", num_reps =", num_reps, "\n")
cat("Grid size:", grid_size, "x", grid_size, "\n")
cat("Starting...\n\n")

# Set seed for reproducibility

set.seed(42)
# Time the original version
time_original = system.time({
  longitudinal_df_keep = c()
  source("MISC/RUN_REPS.R")
  results_original = longitudinal_df_keep
})



cat("\n✓ ORIGINAL VERSION COMPLETED\n")
cat("Total time:", round(time_original["elapsed"], 2), "seconds\n")
cat("User time: ", round(time_original["user.self"], 2), "s\n")
cat("System time:", round(time_original["sys.self"], 2), "s\n")

# ============================================================================

# RUN OPTIMIZED VERSION WITH TIMING

# ============================================================================

cat("\n")
cat("=" = rep("=", 70), sep="")
cat("\nRUNNING OPTIMIZED VERSION (RUN_REPS_OPTIMIZED.R)\n")
cat("=" = rep("=", 70), sep="")
cat("\n")
cat("Parameters: t_max =", t_max, ", num_reps =", num_reps, "\n")
cat("Grid size:", grid_size, "x", grid_size, "\n")
cat("Starting...\n\n")



# CRITICAL: Reset seed to same value for fair comparison

set.seed(42)
# Time the optimized version

time_optimized = system.time({
  longitudinal_df_keep = c()
  source("MISC/RUN_REPS_OPTIMIZED.R")
  results_optimized = longitudinal_df_keep
  
})


cat("\n✓ OPTIMIZED VERSION COMPLETED\n")
cat("Total time:", round(time_optimized["elapsed"], 2), "seconds\n")
cat("User time: ", round(time_optimized["user.self"], 2), "s\n")
cat("System time:", round(time_optimized["sys.self"], 2), "s\n")

# ============================================================================
# PERFORMANCE COMPARISON
# ============================================================================
cat("\n")
cat("=" = rep("=", 70), sep="")
cat("\n         PERFORMANCE COMPARISON RESULTS\n")
cat("=" = rep("=", 70), sep="")
cat("\n\n")

elapsed_original  = time_original["elapsed"]
elapsed_optimized = time_optimized["elapsed"]
speedup = elapsed_original / elapsed_optimized
time_saved = elapsed_original - elapsed_optimized
percent_faster = ((elapsed_original - elapsed_optimized) / elapsed_original) * 100

cat("ORIGINAL VERSION:  ", sprintf("%6.2f", elapsed_original), "seconds\n")
cat("OPTIMIZED VERSION: ", sprintf("%6.2f", elapsed_optimized), "seconds\n")
cat("\n")
cat("SPEEDUP:           ", sprintf("%6.2f", speedup), "x faster\n")
cat("TIME SAVED:        ", sprintf("%6.2f", time_saved), "seconds (",
    sprintf("%.1f", percent_faster), "% faster)\n")

if (speedup >= 3) {
  cat("\n🚀 EXCELLENT! Optimization provides significant speedup!\n")
} else if (speedup >= 1.5) {
  cat("\n✓ GOOD! Noticeable performance improvement.\n")
} else if (speedup >= 1.1) {
  cat("\n✓ Modest improvement. May be more significant with larger t_max.\n")
} else {
  cat("\n⚠ Limited improvement. Try increasing t_max for better assessment.\n")
}

