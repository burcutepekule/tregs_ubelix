# Test to verify update_ROS_batch_cpp behavior
library(Rcpp)

# Source the C++ functions
source("./MISC/FAST_FUNCTIONS_CPP.R")

# Create test data
set.seed(123)
grid_size <- 10
n_phagocytes <- 5

# Initialize
ROS_test <- matrix(0, grid_size, grid_size)
phagocyte_x <- c(2, 5, 7, 2, 9)  # Note: phagocytes 1 and 4 both at x=2
phagocyte_y <- c(3, 6, 8, 3, 4)  # Note: phagocytes 1 and 4 both at y=3
phagocyte_phenotype <- c(1, 0, 1, 1, 0)  # Phagocytes 1, 3, 4 are M1
phagocyte_activity_ROS <- c(0.5, 0.2, 0.3, 0.4, 0.1)
add_ROS <- 0.1

M1_phagocytes <- which(phagocyte_phenotype == 1)
cat("M1 phagocytes:", M1_phagocytes, "\n")
cat("Their positions:\n")
for (i in M1_phagocytes) {
  cat("  Phagocyte", i, "at (", phagocyte_x[i], ",", phagocyte_y[i], ") with activity",
      phagocyte_activity_ROS[i], "\n")
}

# Test C++ version
if (exists("update_ROS_batch_cpp", mode = "function")) {
  ROS_cpp <- update_ROS_batch_cpp(
    ROS_test, M1_phagocytes, phagocyte_x, phagocyte_y,
    phagocyte_activity_ROS, add_ROS
  )
  cat("\nC++ version result:\n")
  cat("  ROS[3,2] =", ROS_cpp[3, 2], "(should be", (0.5 + 0.4) * add_ROS, "from 2 phagocytes)\n")
  cat("  ROS[8,7] =", ROS_cpp[8, 7], "(should be", 0.3 * add_ROS, "from 1 phagocyte)\n")
  cat("  Total ROS =", sum(ROS_cpp), "\n")
} else {
  cat("\nC++ function not available!\n")
}

# Test R version
ROS_R <- matrix(0, grid_size, grid_size)
for (i in M1_phagocytes) {
  ROS_R[phagocyte_y[i], phagocyte_x[i]] <- ROS_R[phagocyte_y[i], phagocyte_x[i]] +
    phagocyte_activity_ROS[i] * add_ROS
}
cat("\nR version result:\n")
cat("  ROS[3,2] =", ROS_R[3, 2], "(should be", (0.5 + 0.4) * add_ROS, "from 2 phagocytes)\n")
cat("  ROS[8,7] =", ROS_R[8, 7], "(should be", 0.3 * add_ROS, "from 1 phagocyte)\n")
cat("  Total ROS =", sum(ROS_R), "\n")

# Compare
if (exists("update_ROS_batch_cpp", mode = "function")) {
  cat("\n=== COMPARISON ===\n")
  if (identical(ROS_cpp, ROS_R)) {
    cat("✓ Results MATCH perfectly!\n")
  } else {
    cat("✗ Results DIFFER!\n")
    diff_matrix <- ROS_cpp - ROS_R
    cat("Max difference:", max(abs(diff_matrix)), "\n")
    cat("Positions with differences:\n")
    diff_idx <- which(abs(diff_matrix) > 1e-10, arr.ind = TRUE)
    if (nrow(diff_idx) > 0) {
      for (i in 1:nrow(diff_idx)) {
        r <- diff_idx[i, 1]
        c <- diff_idx[i, 2]
        cat("  [", r, ",", c, "] C++:", ROS_cpp[r, c], "R:", ROS_R[r, c],
            "Diff:", diff_matrix[r, c], "\n")
      }
    }
  }
}
