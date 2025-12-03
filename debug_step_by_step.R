# ============================================================================
# STEP-BY-STEP DIAGNOSTIC: Track where C++ and R diverge
# ============================================================================
# Run the first few timesteps with both versions and compare state
# ============================================================================

library(dplyr)

# Use the exact same parameters as your main script
grid_size <- 25
n_phagocytes <- round(grid_size * grid_size * 0.05)
n_tregs <- round(grid_size * grid_size * 0.05)
n_commensals_lp <- 20
injury_percentage <- 60
max_level_injury <- 5

# Simple fixed parameters
k_in <- 0.044
x0_in <- 50
act_radius_ROS <- 1
act_radius_DAMPs <- 1
act_radius_SAMPs <- 1

# Load parameters from CSV
params_df <- read.csv("./lhs_parameters_ubelix_merged.csv", stringsAsFactors = FALSE)
param_set_use <- params_df %>% dplyr::filter(param_set_id == 5)

# Source parameter assignment
sterile <- 1
allow_tregs <- 1
randomize_tregs <- 0
source("./MISC/ASSIGN_PARAMETERS.R")

# Function to capture state
capture_state <- function(t, epithelium, phagocytes_df, tregs_df,
                         pathogen_coords, commensal_coords,
                         DAMPs, SAMPs, ROS) {
  list(
    t = t,
    epithelium_injury_sum = sum(epithelium$level_injury),
    epithelium_injury_max = max(epithelium$level_injury),
    n_M0 = sum(phagocytes_df$phenotype == 0),
    n_M1 = sum(phagocytes_df$phenotype == 1),
    n_M2 = sum(phagocytes_df$phenotype == 2),
    n_tregs_active = sum(tregs_df$phenotype == 1),
    n_pathogens = nrow(pathogen_coords),
    n_commensals = nrow(commensal_coords),
    DAMPs_sum = sum(DAMPs),
    SAMPs_sum = sum(SAMPs),
    ROS_sum = sum(ROS),
    DAMPs_max = max(DAMPs),
    SAMPs_max = max(SAMPs),
    ROS_max = max(ROS)
  )
}

# Function to compare states
compare_states <- function(state_R, state_cpp, tolerance = 1e-10) {
  diffs <- list()
  for (name in names(state_R)) {
    if (name == "t") next
    diff <- abs(state_R[[name]] - state_cpp[[name]])
    if (diff > tolerance) {
      diffs[[name]] <- list(
        R = state_R[[name]],
        cpp = state_cpp[[name]],
        diff = diff
      )
    }
  }
  diffs
}

# ============================================================================
# RUN WITH R VERSION (cpp_on = F)
# ============================================================================
cat("\n=== Running with R version ===\n")
set.seed(12345)
source("./MISC/FAST_FUNCTIONS.R")

# Initialize (copy from RUN_REPS_CPP.R)
DAMPs_R <- matrix(0, grid_size, grid_size)
SAMPs_R <- matrix(0, grid_size, grid_size)
ROS_R <- matrix(0, grid_size, grid_size)

epithelium_R <- data.frame(
  x = seq(1, grid_size, 1),
  y = rep(0, grid_size),
  level_injury = 0,
  id = seq(1, grid_size)
)
epithelium_R[injury_site, ]$level_injury <- 1

# Store initial state
states_R <- list()
states_R[[1]] <- list(
  t = 0,
  epithelium_injury_sum = sum(epithelium_R$level_injury),
  DAMPs_sum = sum(DAMPs_R),
  SAMPs_sum = sum(SAMPs_R),
  ROS_sum = sum(ROS_R)
)

cat("R version initialized:\n")
cat("  Injury sum:", states_R[[1]]$epithelium_injury_sum, "\n")
cat("  Random state:", paste(.Random.seed[1:5], collapse=","), "\n\n")

# ============================================================================
# RUN WITH C++ VERSION (cpp_on = T)
# ============================================================================
cat("\n=== Running with C++ version ===\n")
set.seed(12345)
source("./MISC/FAST_FUNCTIONS_CPP.R")

# Initialize (same as R)
DAMPs_cpp <- matrix(0, grid_size, grid_size)
SAMPs_cpp <- matrix(0, grid_size, grid_size)
ROS_cpp <- matrix(0, grid_size, grid_size)

epithelium_cpp <- data.frame(
  x = seq(1, grid_size, 1),
  y = rep(0, grid_size),
  level_injury = 0,
  id = seq(1, grid_size)
)
epithelium_cpp[injury_site, ]$level_injury <- 1

states_cpp <- list()
states_cpp[[1]] <- list(
  t = 0,
  epithelium_injury_sum = sum(epithelium_cpp$level_injury),
  DAMPs_sum = sum(DAMPs_cpp),
  SAMPs_sum = sum(SAMPs_cpp),
  ROS_sum = sum(ROS_cpp)
)

cat("C++ version initialized:\n")
cat("  Injury sum:", states_cpp[[1]]$epithelium_injury_sum, "\n")
cat("  Random state:", paste(.Random.seed[1:5], collapse=","), "\n\n")

# ============================================================================
# COMPARE
# ============================================================================
cat("\n=== Initial State Comparison ===\n")
if (identical(states_R[[1]], states_cpp[[1]])) {
  cat("✓ Initial states match\n")
} else {
  cat("✗ Initial states DIFFER!\n")
  print("R:")
  print(states_R[[1]])
  print("C++:")
  print(states_cpp[[1]])
}

cat("\n=== DIAGNOSTIC COMPLETE ===\n")
cat("Now run a few timesteps of your actual simulation with both cpp_on=T and cpp_on=F\n")
cat("and compare the outputs at each timestep to find where they diverge.\n")
