rm(list=ls())
library(dplyr)
library(tidyr)
library(lhs)
library(readr)

# Define the original parameter bounds (matching mass_simulation_LHS.py)
param_bounds = list(
  th_ROS_microbe = c(0, 1),
  th_ROS_epith_recover = c(0, 1),
  epith_recovery_chance = c(0, 1),
  rat_com_pat_threshold = c(0.5, 1),
  diffusion_speed_DAMPs = c(0, 0.12),
  diffusion_speed_SAMPs = c(0, 0.12),
  diffusion_speed_ROS = c(0, 0.12),
  add_ROS = c(0, 1),
  add_DAMPs = c(0, 1),
  add_SAMPs = c(0, 1),
  ros_decay = c(0, 1),
  DAMPs_decay = c(0, 1),
  SAMPs_decay = c(0, 1),
  activation_threshold_DAMPs = c(0, 1),
  activation_threshold_SAMPs = c(0, 1),
  activity_engulf_M0_baseline = c(0, 0.5),
  activity_engulf_M1_baseline = c(0, 0.5),
  activity_engulf_M2_baseline = c(0, 0.5),
  activity_ROS_M1_baseline = c(0, 0.5),
  rate_leak_commensal_injury = c(0.5, 1),
  rate_leak_pathogen_injury = c(0.5, 1),
  rate_leak_commensal_baseline = c(0, 0.25),
  active_age_limit = c(3, 30),  # discrete parameter, will be rounded
  treg_discrimination_efficiency = c(0, 1)
)

param_names = names(param_bounds)

# Set parameters
set.seed(123)
n_samples = 1e5  # Total number of LHS samples to generate

# Generate LHS samples from original bounds
cat("Generating", n_samples, "LHS samples from original parameter bounds...\n")

n_params = length(param_names)
lhs_unit = randomLHS(n_samples, n_params)

# Create dataframe for samples
lhs_samples = as.data.frame(lhs_unit)
names(lhs_samples) = param_names

# Scale each parameter to its original bounds
for (param in param_names) {
  param_min = param_bounds[[param]][1]
  param_max = param_bounds[[param]][2]
  lhs_samples[[param]] = lhs_samples[[param]] * (param_max - param_min) + param_min
}

# Round the discrete parameter
lhs_samples$active_age_limit = round(lhs_samples$active_age_limit)

# Add parameter set ID
lhs_samples$param_set_id = 0:(nrow(lhs_samples) - 1)

# Reorder columns to put param_set_id first
lhs_samples = lhs_samples[c('param_set_id', param_names)]

# Export
output_file = "/storage/homefs/bt25p365/tregs/lhs_parameters_ubelix.csv"
write.csv(lhs_samples, output_file, row.names = FALSE)
cat("\nDataset saved to:", output_file, "\n")
cat("Total parameter sets:", nrow(lhs_samples), "\n")
