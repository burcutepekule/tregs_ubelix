# ============================================================================
# DIAGNOSTIC SCRIPT TO COMPARE C++ VS R IMPLEMENTATIONS
# ============================================================================
# This script runs a minimal simulation with both implementations and
# compares outputs at each time step to pinpoint where differences occur
# ============================================================================

rm(list = ls())
library(dplyr)

# ============================================================================
# PARAMETERS
# ============================================================================
set.seed(42)
grid_size <- 10
t_max <- 10  # Short run for debugging
n_phagocytes <- 5
n_tregs <- 5
n_commensals_lp <- 5
n_pathogens_lp <- 2

# Fixed parameters
k_in <- 0.044
x0_in <- 50
max_cell_value_ROS <- 1
max_cell_value_DAMPs <- 1
max_cell_value_SAMPs <- 1
act_radius_ROS <- 1
act_radius_DAMPs <- 1
act_radius_SAMPs <- 1

# ============================================================================
# TEST EACH C++ FUNCTION AGAINST R FALLBACK
# ============================================================================

cat("\n=== TESTING C++ FUNCTIONS ===\n\n")

# Load C++ functions
cpp_on <- TRUE
if (cpp_on) {
  source("./MISC/FAST_FUNCTIONS_CPP.R")
} else {
  source("./MISC/FAST_FUNCTIONS.R")
}

# ============================================================================
# Test 1: diffuse_matrix
# ============================================================================
cat("Test 1: diffuse_matrix\n")
test_mat <- matrix(runif(100), 10, 10)
D <- 0.1
max_val <- 1.0

# R version
source("./MISC/FAST_FUNCTIONS.R")
result_R <- diffuse_matrix(test_mat, D, max_val)

# C++ version (if available)
source("./MISC/FAST_FUNCTIONS_CPP.R")
if (exists("diffuse_matrix_cpp", mode = "function")) {
  result_cpp <- diffuse_matrix_cpp(test_mat, D, max_val)
  diff <- max(abs(result_R - result_cpp))
  cat(sprintf("  Max difference: %.15f\n", diff))
  if (diff > 1e-10) {
    cat("  WARNING: Significant difference detected!\n")
    cat("  R result summary:\n")
    print(summary(as.vector(result_R)))
    cat("  C++ result summary:\n")
    print(summary(as.vector(result_cpp)))
  } else {
    cat("  ✓ Results match\n")
  }
} else {
  cat("  C++ version not available\n")
}

# ============================================================================
# Test 2: get_8n_avg_signal
# ============================================================================
cat("\nTest 2: get_8n_avg_signal_vectorized\n")
signal_matrix <- matrix(runif(100), 10, 10)
x_vec <- c(3, 5, 7)
y_vec <- c(2, 5, 8)
act_radius <- 1

# R version
source("./MISC/FAST_FUNCTIONS.R")
result_R <- get_8n_avg_signal_vectorized(x_vec, y_vec, act_radius, signal_matrix, grid_size)

# C++ version
source("./MISC/FAST_FUNCTIONS_CPP.R")
if (exists("get_8n_avg_signal_vectorized_cpp", mode = "function")) {
  result_cpp <- get_8n_avg_signal_vectorized_cpp(x_vec, y_vec, act_radius, signal_matrix, grid_size)
  diff <- max(abs(result_R - result_cpp))
  cat(sprintf("  Max difference: %.15f\n", diff))
  if (diff > 1e-10) {
    cat("  WARNING: Significant difference detected!\n")
    cat("  R results:", result_R, "\n")
    cat("  C++ results:", result_cpp, "\n")
  } else {
    cat("  ✓ Results match\n")
  }
} else {
  cat("  C++ version not available\n")
}

# ============================================================================
# Test 3: iszero_coordinates (R only - checking for wrapper)
# ============================================================================
cat("\nTest 3: iszero_coordinates\n")
set.seed(123)
x_test <- c(0, 1, -1, 0, 1)

# R version
source("./MISC/FAST_FUNCTIONS.R")
set.seed(123)
result_R <- iszero_coordinates(x_test)

# Check if there's a C++ wrapper being called
source("./MISC/FAST_FUNCTIONS_CPP.R")
set.seed(123)
result_with_cpp_loaded <- iszero_coordinates(x_test)

if (!identical(result_R, result_with_cpp_loaded)) {
  cat("  WARNING: iszero_coordinates gives different results when C++ is loaded!\n")
  cat("  R only:", result_R, "\n")
  cat("  With C++ loaded:", result_with_cpp_loaded, "\n")
} else {
  cat("  ✓ Same function being used (no C++ wrapper)\n")
}

# ============================================================================
# Test 4: shift_insert_fast
# ============================================================================
cat("\nTest 4: shift_insert_fast\n")
vec <- c(1, 2, 3, 4, 5)
insert_vals <- c(9, 8)

# R version
source("./MISC/FAST_FUNCTIONS.R")
result_R <- shift_insert_fast(vec, insert_vals)

# C++ version
source("./MISC/FAST_FUNCTIONS_CPP.R")
if (exists("shift_insert_fast_cpp", mode = "function")) {
  result_cpp <- shift_insert_fast_cpp(vec, insert_vals)
  if (!identical(result_R, result_cpp)) {
    cat("  WARNING: Results differ!\n")
    cat("  R result:", result_R, "\n")
    cat("  C++ result:", result_cpp, "\n")
  } else {
    cat("  ✓ Results match\n")
  }
} else {
  cat("  C++ version not available\n")
}

# ============================================================================
# Test 5: Check which C++ functions are actually loaded
# ============================================================================
cat("\n=== C++ FUNCTION AVAILABILITY ===\n")
check_cpp_status()

cat("\n=== DIAGNOSIS COMPLETE ===\n")
cat("If any differences were found above, that's likely the source of the discrepancy.\n")
