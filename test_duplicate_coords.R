# Test: Do R and C++ handle multiple agents at same location differently?

# Create a test matrix
test_mat <- matrix(0, 5, 5)

# Simulate 3 tregs at positions, with 2 at the same location
treg_x <- c(2, 3, 2)  # Note: tregs 1 and 3 are both at x=2
treg_y <- c(2, 3, 2)  # Note: tregs 1 and 3 are both at y=2
treg_activity <- c(0.5, 0.5, 0.5)
active_tregs <- c(1, 2, 3)

cat("Testing duplicate coordinate handling...\n")
cat("Treg locations: (2,2), (3,3), (2,2)\n")
cat("All tregs contribute 0.5\n\n")

# ============================================================================
# R version (as used in RUN_REPS_CPP.R fallback)
# ============================================================================
test_mat_R <- test_mat
coords = cbind(treg_y[active_tregs], treg_x[active_tregs])
cat("R version coords matrix:\n")
print(coords)
test_mat_R[coords] = test_mat_R[coords] + treg_activity[active_tregs]

cat("\nR version result at (2,2):", test_mat_R[2,2], "\n")
cat("Expected if both tregs add: 1.0\n")
cat("Expected if only last adds: 0.5\n")

# ============================================================================
# C++ version behavior (simulate what C++ does)
# ============================================================================
test_mat_cpp <- test_mat
for (i in 1:length(active_tregs)) {
  idx <- active_tregs[i]
  x <- treg_x[idx]
  y <- treg_y[idx]
  test_mat_cpp[y, x] <- test_mat_cpp[y, x] + treg_activity[idx]
}

cat("\nC++ version result at (2,2):", test_mat_cpp[2,2], "\n")

# ============================================================================
# Compare
# ============================================================================
cat("\n=== RESULT ===\n")
if (test_mat_R[2,2] != test_mat_cpp[2,2]) {
  cat("❌ DIFFERENCE FOUND!\n")
  cat("R version: ", test_mat_R[2,2], "\n")
  cat("C++ version: ", test_mat_cpp[2,2], "\n")
  cat("\nThis explains the cpp_on=T vs cpp_on=F discrepancy!\n")
  cat("When multiple agents are at the same location,\n")
  cat("R matrix indexing only applies the LAST update,\n")
  cat("while C++ correctly accumulates ALL updates.\n")
} else {
  cat("✓ Both versions give same result\n")
  cat("This is not the source of the discrepancy\n")
}
