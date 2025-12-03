# Create the complete dataframe with bacteria registry
phagocytes <- data.frame(
  x = phagocyte_x,
  y = phagocyte_y,
  pathogens_engulfed = phagocyte_pathogens_engulfed,
  commensals_engulfed = phagocyte_commensals_engulfed,
  num_times_activated = phagocyte_num_times_activated,
  phenotype = phagocyte_phenotype,
  activity_ROS = phagocyte_activity_ROS,
  activity_engulf = phagocyte_activity_engulf,
  active_age = phagocyte_active_age
)

# Add the bacteria registry as a list column
phagocytes$bacteria_registry <- 1

tregs=data.frame(
  x = treg_x,
  y = treg_y,
  active_age = treg_active_age,
  phenotype = treg_phenotype,
  activity_SAMPs_binary = treg_activity_SAMPs_binary
)

if(nrow(pathogen_coords) == 0) {
  pathogens_lp = data.frame(x = numeric(0), y = numeric(0), id = numeric(0))
}else{
  pathogens_lp =   data.frame(
    x = pathogen_coords[, "x"],
    y = pathogen_coords[, "y"]
  )
}

if(nrow(commensal_coords) == 0) {
  commensals_lp = data.frame(x = numeric(0), y = numeric(0), id = numeric(0))
}else{
  commensals_lp =   data.frame(
    x = commensal_coords[, "x"],
    y = commensal_coords[, "y"]
  )
}
