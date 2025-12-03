# ============================================================================
# TARGETED DIAGNOSTIC: Find where C++ and R implementations diverge
# ============================================================================
# Run minimal simulations and compare state at each timestep
# ============================================================================

rm(list = ls())

# Set the same seed for both
SEED <- 12345

# Minimal parameters for debugging
grid_size <- 10
t_max <- 5  # Very short
n_phagocytes <- 3
n_tregs <- 3
n_commensals_lp <- 3

# Fixed parameters (minimal)
k_in <- 0.044
x0_in <- 50
injury_percentage <- 60
max_level_injury <- 5
act_radius_ROS <- 1
act_radius_DAMPs <- 1
act_radius_SAMPs <- 1
act_radius_treg <- 1

# ============================================================================
# Run with R version
# ============================================================================
cat("\n=== Running with R version (cpp_on = F) ===\n")
set.seed(SEED)
source("./MISC/FAST_FUNCTIONS.R")

# Simple test: create a signal matrix and test get_8n_avg_signal
test_matrix <- matrix(runif(100), 10, 10)

# Test at boundary
x_test <- 1
y_test <- 1
result_R <- get_8n_avg_signal_fast(x_test, y_test, act_radius_ROS, test_matrix)
cat("R get_8n_avg_signal at (1,1):", result_R, "\n")

# Test vectorized version
x_vec <- c(1, 5, 10)
y_vec <- c(1, 5, 10)
results_R_vec <- get_8n_avg_signal_vectorized(x_vec, y_vec, act_radius_ROS, test_matrix, grid_size)
cat("R vectorized at (1,1), (5,5), (10,10):", results_R_vec, "\n")

# ============================================================================
# Run with C++ version
# ============================================================================
cat("\n=== Running with C++ version (cpp_on = T) ===\n")
set.seed(SEED)
source("./MISC/FAST_FUNCTIONS_CPP.R")

# Test the same
result_cpp <- if (exists("get_8n_avg_signal_cpp", mode = "function")) {
  get_8n_avg_signal_cpp(x_test, y_test, act_radius_ROS, test_matrix, grid_size)
} else {
  cat("C++ function not available!\n")
  get_8n_avg_signal_fast(x_test, y_test, act_radius_ROS, test_matrix)
}
cat("C++ get_8n_avg_signal at (1,1):", result_cpp, "\n")

results_cpp_vec <- if (exists("get_8n_avg_signal_vectorized_cpp", mode = "function")) {
  get_8n_avg_signal_vectorized_cpp(x_vec, y_vec, act_radius_ROS, test_matrix, grid_size)
} else {
  cat("C++ vectorized function not available!\n")
  get_8n_avg_signal_vectorized(x_vec, y_vec, act_radius_ROS, test_matrix, grid_size)
}
cat("C++ vectorized at (1,1), (5,5), (10,10):", results_cpp_vec, "\n")

# ============================================================================
# Compare
# ============================================================================
cat("\n=== Comparison ===\n")
cat("Difference at (1,1):", abs(result_R - result_cpp), "\n")
cat("Max difference vectorized:", max(abs(results_R_vec - results_cpp_vec)), "\n")

if (max(abs(results_R_vec - results_cpp_vec)) > 1e-10) {
  cat("\nWARNING: Significant difference detected in signal calculation!\n")
  cat("R results:\n")
  print(results_R_vec)
  cat("C++ results:\n")
  print(results_cpp_vec)
  cat("Differences:\n")
  print(results_R_vec - results_cpp_vec)
} else {
  cat("\n✓ Signal calculations match between R and C++\n")
}

# ============================================================================
# Test diffusion
# ============================================================================
cat("\n=== Testing diffusion ===\n")
set.seed(SEED)
test_mat <- matrix(runif(100), 10, 10)
D <- 0.1
max_val <- 1.0

source("./MISC/FAST_FUNCTIONS.R")
result_R_diff <- diffuse_matrix(test_mat, D, max_val)

source("./MISC/FAST_FUNCTIONS_CPP.R")
result_cpp_diff <- if (exists("diffuse_matrix_cpp", mode = "function")) {
  diffuse_matrix_cpp(test_mat, D, max_val)
} else {
  diffuse_matrix(test_mat, D, max_val)
}

diff_max <- max(abs(result_R_diff - result_cpp_diff))
cat("Max difference in diffusion:", sprintf("%.15f", diff_max), "\n")

if (diff_max > 1e-10) {
  cat("WARNING: Diffusion results differ!\n")
} else {
  cat("✓ Diffusion matches\n")
}

# ============================================================================
# Print summary
# ============================================================================
cat("\n=== SUMMARY ===\n")
cat("If all tests passed, the basic C++ functions are working correctly.\n")
cat("If differences persist in your full simulation, the issue may be:\n")
cat("1. Different random number consumption patterns\n")
cat("2. Edge cases not covered by these simple tests\n")
cat("3. Interaction effects between multiple function calls\n")
cat("4. Differences in how specific parameters are handled\n")
