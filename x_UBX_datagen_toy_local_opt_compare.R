rm(list=ls())
library(dplyr)
library(tidyr)

# ============================================================================
# TIMING COMPARISON: ORIGINAL vs OPTIMIZED RUN_REPS.R
# ============================================================================
# This script runs both versions with identical parameters and reports timing

# ============================================================================
# LOAD FUNCTIONS
# ============================================================================
source("MISC/FAST_FUNCTIONS.R")
# source("MISC/PLOT_FUNCTIONS.R")  # Comment out if not needed for speed test

# ============================================================================
# FIXED PARAMETERS (toy example for faster testing)
# ============================================================================
t_max      = 50      # Increase this for more realistic timing (100-500)
num_reps   = 5       # Number of replicates per version
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

# ============================================================================
# SET TEST PARAMETERS
# ============================================================================
# Use a single parameter set for controlled comparison
sterile         = 0
allow_tregs     = 1
randomize_tregs = 0

# Load a test parameter set (or define manually)
params_df = read.csv("lhs_parameters_ubelix.csv", stringsAsFactors = FALSE)
param_set_id_use = params_df$param_set_id[1]  # Use first parameter set
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

# ============================================================================
# VERIFY RESULTS ARE IDENTICAL
# ============================================================================
cat("\n")
cat("=" = rep("=", 70), sep="")
cat("\n         VERIFYING RESULTS ARE IDENTICAL\n")
cat("=" = rep("=", 70), sep="")
cat("\n\n")

# Check if results are identical
results_match = all.equal(results_original, results_optimized)

if (isTRUE(results_match)) {
  cat("✓ PERFECT! Results are IDENTICAL between versions.\n")
  cat("  The optimization maintains 100% logical equivalence.\n")
} else {
  cat("⚠ WARNING: Results differ between versions!\n")
  cat("  Differences found:\n")
  print(results_match)
  cat("\n  This may indicate a bug in the optimization.\n")
}

# Additional numerical comparison
max_diff = max(abs(as.matrix(results_original[, -c(1:6)]) -
                   as.matrix(results_optimized[, -c(1:6)])), na.rm = TRUE)
cat("\nMaximum numerical difference:", format(max_diff, scientific = TRUE), "\n")

if (max_diff < 1e-10) {
  cat("✓ Numerical precision is excellent (< 1e-10)\n")
}

# ============================================================================
# EXTRAPOLATED PERFORMANCE FOR FULL SIMULATION
# ============================================================================
cat("\n")
cat("=" = rep("=", 70), sep="")
cat("\n         EXTRAPOLATED PERFORMANCE FOR FULL SIMULATION\n")
cat("=" = rep("=", 70), sep="")
cat("\n\n")

# Typical full simulation parameters
full_t_max = 5000
full_num_reps = 100
full_scenarios = 6  # Number of scenarios in typical run

# Extrapolate times
scale_factor = (full_t_max / t_max) * (full_num_reps / num_reps)
estimated_original_full = elapsed_original * scale_factor * full_scenarios
estimated_optimized_full = elapsed_optimized * scale_factor * full_scenarios
estimated_savings_full = estimated_original_full - estimated_optimized_full

cat("For a full simulation (t_max =", full_t_max, ", num_reps =", full_num_reps,
    ", scenarios =", full_scenarios, "):\n\n")

cat("Estimated ORIGINAL time:  ", sprintf("%8.1f", estimated_original_full),
    " seconds (", sprintf("%.1f", estimated_original_full/3600), " hours)\n")
cat("Estimated OPTIMIZED time: ", sprintf("%8.1f", estimated_optimized_full),
    " seconds (", sprintf("%.1f", estimated_optimized_full/3600), " hours)\n")
cat("\n")
cat("Estimated TIME SAVINGS:   ", sprintf("%8.1f", estimated_savings_full),
    " seconds (", sprintf("%.1f", estimated_savings_full/3600), " hours)\n")

if (estimated_savings_full > 3600) {
  cat("\n🎉 Using the optimized version could save you ",
      sprintf("%.1f", estimated_savings_full/3600), " hours!\n")
  cat("   DEFINITELY WORTH IT for production runs!\n")
} else if (estimated_savings_full > 600) {
  cat("\n✓ Using the optimized version could save you ",
      sprintf("%.1f", estimated_savings_full/60), " minutes!\n")
  cat("   Recommended for production runs.\n")
}

# ============================================================================
# SUMMARY AND RECOMMENDATIONS
# ============================================================================
cat("\n")
cat("=" = rep("=", 70), sep="")
cat("\n         SUMMARY AND RECOMMENDATIONS\n")
cat("=" = rep("=", 70), sep="")
cat("\n\n")

if (speedup >= 2 && isTRUE(results_match)) {
  cat("✓ RECOMMENDATION: USE THE OPTIMIZED VERSION\n")
  cat("  - Provides ", sprintf("%.1f", speedup), "x speedup\n")
  cat("  - Results are identical\n")
  cat("  - No downsides detected\n\n")
  cat("  To use in production, edit UBX_datagen.R line 100:\n")
  cat("  source(\"MISC/RUN_REPS_OPTIMIZED.R\")\n")
} else if (speedup >= 1.2 && isTRUE(results_match)) {
  cat("✓ RECOMMENDATION: Consider using optimized version\n")
  cat("  - Provides ", sprintf("%.1f", speedup), "x speedup\n")
  cat("  - Results are identical\n")
  cat("  - May be more beneficial with larger simulations\n")
} else {
  cat("⚠ RECOMMENDATION: Further testing needed\n")
  if (!isTRUE(results_match)) {
    cat("  - Results differ between versions - investigate!\n")
  }
  if (speedup < 1.2) {
    cat("  - Limited speedup with current parameters\n")
    cat("  - Try testing with larger t_max (e.g., 500-1000)\n")
  }
}

cat("\n")
cat("=" = rep("=", 70), sep="")
cat("\n")
