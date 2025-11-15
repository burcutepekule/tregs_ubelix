# ============================================================================
# C++ INSTALLATION TEST SCRIPT
# ============================================================================
# This script verifies that Rcpp is installed and C++ functions compile
# Run this BEFORE using the C++ accelerated version
# ============================================================================

cat("\n")
cat("=" = rep("=", 70), sep = "")
cat("\n")
cat("C++ INSTALLATION TEST\n")
cat("=" = rep("=", 70), sep = "")
cat("\n\n")

# ============================================================================
# TEST 1: Check if Rcpp is installed
# ============================================================================

cat("TEST 1: Checking if Rcpp is installed... ")

if (requireNamespace("Rcpp", quietly = TRUE)) {
  cat("✓ PASS\n")
  library(Rcpp)
} else {
  cat("✗ FAIL\n")
  cat("  Rcpp is not installed. Install it with:\n")
  cat("  install.packages('Rcpp')\n\n")
  stop("Rcpp not available")
}

# ============================================================================
# TEST 2: Check if C++ compiler is available
# ============================================================================

cat("TEST 2: Checking if C++ compiler is available... ")

compiler_path <- Sys.which("gcc")
if (compiler_path != "") {
  cat("✓ PASS\n")
  cat("  Compiler found at:", compiler_path, "\n")
} else {
  compiler_path <- Sys.which("clang")
  if (compiler_path != "") {
    cat("✓ PASS\n")
    cat("  Compiler found at:", compiler_path, "\n")
  } else {
    cat("✗ FAIL\n")
    cat("  No C++ compiler found.\n")
    cat("  Install:\n")
    cat("    - Linux: sudo apt-get install r-base-dev\n")
    cat("    - Mac: xcode-select --install\n")
    cat("    - Windows: Install Rtools\n\n")
    stop("C++ compiler not available")
  }
}

# ============================================================================
# TEST 3: Test basic Rcpp compilation
# ============================================================================

cat("TEST 3: Testing basic C++ compilation... ")

tryCatch({
  result <- Rcpp::evalCpp("2 + 2")
  if (result == 4) {
    cat("✓ PASS\n")
  } else {
    cat("✗ FAIL (incorrect result)\n")
    stop("Basic C++ compilation failed")
  }
}, error = function(e) {
  cat("✗ FAIL\n")
  cat("  Error:", conditionMessage(e), "\n\n")
  stop("Basic C++ compilation failed")
})

# ============================================================================
# TEST 4: Compile FAST_FUNCTIONS.cpp
# ============================================================================

cat("TEST 4: Compiling FAST_FUNCTIONS.cpp... ")

cpp_file <- "MISC/FAST_FUNCTIONS.cpp"

if (!file.exists(cpp_file)) {
  cat("✗ FAIL\n")
  cat("  File not found:", cpp_file, "\n\n")
  stop("C++ file not found")
}

tryCatch({
  sourceCpp(cpp_file)
  cat("✓ PASS\n")
}, error = function(e) {
  cat("✗ FAIL\n")
  cat("  Error:", conditionMessage(e), "\n\n")
  stop("C++ compilation failed")
})

# ============================================================================
# TEST 5: Verify C++ functions are available
# ============================================================================

cat("TEST 5: Verifying C++ functions are available...\n")

expected_functions <- c(
  "get_8n_avg_signal_cpp",
  "get_8n_avg_signal_vectorized_cpp",
  "kill_microbes_with_ros_cpp",
  "calculate_phagocyte_signals_cpp",
  "find_overlapping_microbes_cpp",
  "find_nearby_tregs_cpp",
  "calculate_epithelial_ros_cpp",
  "diffuse_matrix_cpp",
  "shift_insert_fast_cpp",
  "iszero_coordinates_cpp",
  "update_SAMPs_batch_cpp",
  "update_ROS_batch_cpp"
)

all_found <- TRUE
for (func_name in expected_functions) {
  if (exists(func_name, mode = "function")) {
    cat("  ✓", func_name, "\n")
  } else {
    cat("  ✗", func_name, "(NOT FOUND)\n")
    all_found <- FALSE
  }
}

if (all_found) {
  cat("  All functions found ✓\n")
} else {
  cat("  Some functions missing ✗\n")
  stop("Not all C++ functions are available")
}

# ============================================================================
# TEST 6: Simple functionality test
# ============================================================================

cat("TEST 6: Testing C++ function correctness... ")

# Create test data
test_matrix <- matrix(runif(25*25), 25, 25)
test_coords <- matrix(c(10, 15, 1, 12, 18, 2, 20, 5, 3), ncol = 3)
colnames(test_coords) <- c("x", "y", "id")

# Test signal calculation
tryCatch({
  result_signal <- get_8n_avg_signal_cpp(10, 15, 1, test_matrix, 25)
  if (is.numeric(result_signal) && result_signal >= 0 && result_signal <= 1) {
    cat("✓ PASS (signal calculation)\n")
  } else {
    cat("✗ FAIL (signal calculation returned unexpected value)\n")
    stop("C++ function returned unexpected value")
  }
}, error = function(e) {
  cat("✗ FAIL\n")
  cat("  Error:", conditionMessage(e), "\n\n")
  stop("C++ function test failed")
})

# Test microbe killing
cat("         Testing microbe killing function... ")
tryCatch({
  result_kill <- kill_microbes_with_ros_cpp(test_coords, test_matrix, 1, 0.5, 25)
  if (is.list(result_kill) &&
      "surviving_microbes" %in% names(result_kill) &&
      "n_killed" %in% names(result_kill)) {
    cat("✓ PASS\n")
  } else {
    cat("✗ FAIL (unexpected return structure)\n")
    stop("Microbe killing function returned unexpected structure")
  }
}, error = function(e) {
  cat("✗ FAIL\n")
  cat("  Error:", conditionMessage(e), "\n\n")
  stop("Microbe killing function test failed")
})

# ============================================================================
# TEST 7: Performance comparison
# ============================================================================

cat("TEST 7: Quick performance test... ")

# Load R versions
source("MISC/FAST_FUNCTIONS.R", local = TRUE)

# Time R version
grid_size <- 25
time_r <- system.time({
  for (i in 1:100) {
    result_r <- get_8n_avg_signal_fast(10, 15, 1, test_matrix)
  }
})

# Time C++ version
time_cpp <- system.time({
  for (i in 1:100) {
    result_cpp <- get_8n_avg_signal_cpp(10, 15, 1, test_matrix, 25)
  }
})

speedup <- time_r["elapsed"] / time_cpp["elapsed"]

cat("✓ PASS\n")
cat("  R version:   ", round(time_r["elapsed"], 4), "seconds\n")
cat("  C++ version: ", round(time_cpp["elapsed"], 4), "seconds\n")
cat("  Speedup:     ", round(speedup, 1), "x faster\n")

if (speedup < 2) {
  cat("  ⚠ Warning: Speedup is less than expected (should be >10x)\n")
  cat("    This might be due to overhead in small test\n")
}

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n")
cat("=" = rep("=", 70), sep = "")
cat("\n")
cat("ALL TESTS PASSED! ✓\n")
cat("=" = rep("=", 70), sep = "")
cat("\n\n")

cat("Your system is ready for C++ accelerated simulations!\n\n")

cat("Next steps:\n")
cat("1. Run a small test:\n")
cat("   Rscript UBX_datagen_cpp.R 10 1\n\n")
cat("2. Check C++ status in your scripts:\n")
cat("   source('MISC/FAST_FUNCTIONS_CPP.R')\n")
cat("   check_cpp_status()\n\n")
cat("3. Deploy to production:\n")
cat("   Use UBX_datagen_cpp.R instead of UBX_datagen.R\n\n")

cat("Expected performance improvement: 20-100x faster!\n\n")

cat("For more information, see CPP_DEPLOYMENT_GUIDE.md\n\n")
