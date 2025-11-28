# vector of parameter names
if(M1_M2_diff==1){
  param_names = c(
    "th_ROS_microbe",
    "th_ROS_epith_recover",
    "epith_recovery_chance",
    "rat_com_pat_threshold",
    "diffusion_speed_DAMPs",
    "diffusion_speed_SAMPs",
    "diffusion_speed_ROS",
    "add_ROS",
    "add_DAMPs",
    "add_SAMPs",
    "ros_decay",
    "DAMPs_decay",
    "SAMPs_decay",
    "activation_threshold_DAMPs",
    "activation_threshold_SAMPs",
    "activity_engulf_M0_baseline",
    # "activity_engulf_M1_baseline", #remove for colinearity
    # "activity_engulf_M2_baseline", #remove for colinearity
    "activity_ROS_M1_baseline",
    "rate_leak_commensal_injury",
    "rate_leak_pathogen_injury",
    "rate_leak_commensal_baseline",
    "active_age_limit",
    "treg_discrimination_efficiency",
    "activity_engulf_M1_M2_diff"
  )
}else{
  param_names = c(
    "th_ROS_microbe",
    "th_ROS_epith_recover",
    "epith_recovery_chance",
    "rat_com_pat_threshold",
    "diffusion_speed_DAMPs",
    "diffusion_speed_SAMPs",
    "diffusion_speed_ROS",
    "add_ROS",
    "add_DAMPs",
    "add_SAMPs",
    "ros_decay",
    "DAMPs_decay",
    "SAMPs_decay",
    "activation_threshold_DAMPs",
    "activation_threshold_SAMPs",
    "activity_engulf_M0_baseline",
    "activity_engulf_M1_baseline",
    "activity_engulf_M2_baseline",
    "activity_ROS_M1_baseline",
    "rate_leak_commensal_injury",
    "rate_leak_pathogen_injury",
    "rate_leak_commensal_baseline",
    "active_age_limit",
    "treg_discrimination_efficiency"
  )
}
