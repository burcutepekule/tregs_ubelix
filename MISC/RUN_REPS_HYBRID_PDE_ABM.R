# ============================================================================
# HYBRID PDE-ABM SIMULATION LOOP
# ============================================================================
# PDE: Microbes (P, C) and signals (DAMPs, SAMPs, ROS) as concentration fields
# ABM: Epithelial cells, macrophages, Tregs as discrete agents
# Memory: Exponential decay for macrophage engulfment history
# ============================================================================

for (reps_in in 0:(num_reps-1)){

  # ==========================================================================
  # INITIALIZE PDE FIELDS (concentration grids)
  # ==========================================================================
  DAMPs = matrix(0, grid_size, grid_size)
  SAMPs = matrix(0, grid_size, grid_size)
  ROS = matrix(0, grid_size, grid_size)

  # Microbe concentrations (NEW: continuous fields instead of discrete agents)
  P_field = matrix(0, grid_size, grid_size)  # Pathogen concentration
  C_field = matrix(0, grid_size, grid_size)  # Commensal concentration

  # ==========================================================================
  # INITIALIZE ABM AGENTS
  # ==========================================================================

  # Epithelial cells (discrete, at y=1)
  epithelium = data.frame(
    x = seq(1, grid_size, 1),
    y = rep(1, grid_size),
    level_injury = 0,
    id = seq(1, grid_size)
  )
  epithelium[injury_site, ]$level_injury = 1

  # Macrophages (discrete agents)
  phagocyte_x = sample(1:grid_size, n_phagocytes, TRUE)
  phagocyte_y = sample(2:grid_size, n_phagocytes, TRUE)
  phagocyte_phenotype = rep(0, n_phagocytes)  # 0=M0, 1=M1, 2=M2
  phagocyte_activity_ROS = rep(activity_ROS_M0, n_phagocytes) # All start from M0 phenotype
  phagocyte_activity_engulf = rep(activity_engulf_M0, n_phagocytes) # All start from M0 phenotype
  phagocyte_active_age = rep(0, n_phagocytes)

  # Decaying memory 
  phagocyte_pathogen_memory = rep(0, n_phagocytes)
  phagocyte_commensal_memory = rep(0, n_phagocytes)
  memory_decay_rate = exp(-log(2) / memory_half_life)

  # Tregs (discrete agents)
  treg_x = sample(1:grid_size, n_tregs, TRUE)
  treg_y = sample(2:grid_size, n_tregs, TRUE)
  treg_active_age = rep(0, n_tregs)
  treg_phenotype = rep(0, n_tregs)  # 0=resting, 1=activated
  treg_activity_SAMPs_binary = rep(0, n_tregs)

  # Initialize microbe fields with starting populations
  # Pathogens at injury sites
  n_pathogens_lp = round(rate_leak_pathogen_injury * length(injury_site))
  if (n_pathogens_lp > 0) {
    for (i in 1:n_pathogens_lp) {
      x_start = sample(injury_site, 1)
      P_field[1, x_start] = P_field[1, x_start] + 1
    }
  }

  # Commensals distributed
  for (i in 1:n_commensals_lp) {
    x_start = sample(1:grid_size, 1)
    y_start = sample(1:grid_size, 1)
    C_field[y_start, x_start] = C_field[y_start, x_start] + 1
  }

  # Death counters (for PDE: track total removal, not individual kills)
  pathogens_killed_by_ROS = 0
  pathogens_killed_by_Mac = rep(0, 3)
  commensals_killed_by_ROS = 0
  commensals_killed_by_Mac = rep(0, 3)

  # Longitudinal tracking
  epithelium_longitudinal = matrix(0, nrow = t_max, ncol = (max_level_injury + 1))
  macrophages_longitudinal = matrix(0, nrow = t_max, ncol = 3)
  microbes_longitudinal = matrix(0, nrow = t_max, ncol = 2)
  tregs_longitudinal = matrix(0, nrow = t_max, ncol = 2)
  microbes_cumdeath_longitudinal = matrix(0, nrow = t_max, ncol = 2 * 4)

  # ==========================================================================
  # MAIN SIMULATION LOOP
  # ==========================================================================
  for (t in 1:t_max) {

    # Update injury site
    injury_site_updated = which(epithelium$level_injury > 0)

    # ========================================================================
    # STEP 1: AGENT ACTIONS → PDE SOURCES
    # ========================================================================

    # 1.1 SAMPs production (from activated Tregs)
    active_tregs = which(treg_phenotype == 1)
    if (length(active_tregs) > 0) {
      SAMPs = update_SAMPs_batch_cpp(
        SAMPs, active_tregs, treg_x, treg_y,
        treg_activity_SAMPs_binary, add_SAMPs, allow_tregs
      )
    }

    # 1.2 ROS production (from M1 macrophages)
    M1_phagocytes = which(phagocyte_phenotype == 1)
    if (length(M1_phagocytes) > 0) {
      coords = cbind(phagocyte_y[M1_phagocytes], phagocyte_x[M1_phagocytes])
      ROS[coords] = ROS[coords] + phagocyte_activity_ROS[M1_phagocytes] * add_ROS
    }

    # 1.3 DAMP production (from epithelial injury)
    DAMPs[1, ] = DAMPs[1, ] + epithelium$level_injury * add_DAMPs

    # 1.4 DAMP production (from microbe contact at epithelium)
    microbe_load_epithelium = P_field[1, ] + C_field[1, ]
    DAMPs[1, ] = DAMPs[1, ] + logistic_scaled_0_to_5_quantized(microbe_load_epithelium) * add_DAMPs

    # ========================================================================
    # STEP 2: PDE EVOLUTION
    # ========================================================================

    # 2.1 Microbe entry from epithelium
    for (k in 1:nrow(epithelium)) {
      x = epithelium$x[k]
      inj = epithelium$level_injury[k]

      # Pathogen leakage
      n_new_pathogens = rate_leak_pathogen_injury*inj
      P_field[1, x] = P_field[1, x] + n_new_pathogens

      # Commensal leakage
      n_new_commensals = rate_leak_commensal_baseline + rate_leak_commensal_injury * inj
      C_field[1, x] = C_field[1, x] + n_new_commensals
    }

    # 2.2 Diffusion for all fields
    P_field = diffuse_matrix_cpp(P_field, D_microbe, max_microbe_concentration)
    C_field = diffuse_matrix_cpp(C_field, D_microbe, max_microbe_concentration)
    DAMPs = diffuse_matrix_cpp(DAMPs, diffusion_speed_DAMPs, max_cell_value_DAMPs)
    SAMPs = diffuse_matrix_cpp(SAMPs, diffusion_speed_SAMPs, max_cell_value_SAMPs)
    ROS = diffuse_matrix_cpp(ROS, diffusion_speed_ROS, max_cell_value_ROS)

    # 2.3 Decay
    DAMPs = DAMPs - DAMPs_decay * DAMPs
    SAMPs = SAMPs - SAMPs_decay * SAMPs
    ROS = ROS - ros_decay * ROS

    # ========================================================================
    # STEP 3: ABM MOVEMENT (Chemotaxis)
    # ========================================================================

    # 3.1 Tregs movement
    density_matrix_tregs = if (randomize_tregs == 1) matrix(1, grid_size, grid_size) else DAMPs
    all_equal_treg = all(density_matrix_tregs == density_matrix_tregs[1, 1])

    if (!all_equal_treg) {
      for (i in 1:length(treg_x)) {
        x = treg_x[i]
        y = treg_y[i]
        x_range = max(1, x - 1):min(grid_size, x + 1)
        y_range = max(1, y - 1):min(grid_size, y + 1)
        neighbors_x = rep(x_range, each = length(y_range))
        neighbors_y = rep(y_range, times = length(x_range))
        neighbor_densities = density_matrix_tregs[cbind(neighbors_y, neighbors_x)]
        total = sum(neighbor_densities)
        if (total > 0) {
          probs = neighbor_densities / total
        } else {
          probs = rep(1 / length(neighbor_densities), length(neighbor_densities))
        }
        chosen_idx = sample(1:length(neighbors_x), 1, prob = probs)
        treg_x[i] = neighbors_x[chosen_idx]
        treg_y[i] = neighbors_y[chosen_idx]
      }
    } else {
      dy_treg = ifelse(treg_y == 1,
                       sample(c(1), size = length(treg_y), replace = TRUE),
                       sample(c(-1, 0, 1), size = length(treg_y), replace = TRUE))
      dx_treg = iszero_coordinates(dy_treg)
      treg_x = pmin(pmax(treg_x + dx_treg, 1), grid_size)
      treg_y = pmin(pmax(treg_y + dy_treg, 1), grid_size)
    }

    # 3.2 Macrophage movement
    density_matrix_phagocytes = DAMPs
    all_equal_phagocytes = all(density_matrix_phagocytes == density_matrix_phagocytes[1, 1])

    if (!all_equal_phagocytes) {
      for (i in 1:length(phagocyte_x)) {
        x = phagocyte_x[i]
        y = phagocyte_y[i]
        x_range = max(1, x - 1):min(grid_size, x + 1)
        y_range = max(1, y - 1):min(grid_size, y + 1)
        neighbors_x = rep(x_range, each = length(y_range))
        neighbors_y = rep(y_range, times = length(x_range))
        neighbor_densities = density_matrix_phagocytes[cbind(neighbors_y, neighbors_x)]
        total = sum(neighbor_densities)
        if (total > 0) {
          probs = neighbor_densities / total
        } else {
          probs = rep(1 / length(neighbor_densities), length(neighbor_densities))
        }
        chosen_idx = sample(1:length(neighbors_x), 1, prob = probs)
        phagocyte_x[i] = neighbors_x[chosen_idx]
        phagocyte_y[i] = neighbors_y[chosen_idx]
      }
    } else {
      dy_phagocyte = ifelse(phagocyte_y == 1,
                            sample(c(1), size = length(phagocyte_y), replace = TRUE),
                            sample(c(-1, 0, 1), size = length(phagocyte_y), replace = TRUE))
      dx_phagocyte = iszero_coordinates(dy_phagocyte)
      phagocyte_x = pmin(pmax(phagocyte_x + dx_phagocyte, 1), grid_size)
      phagocyte_y = pmin(pmax(phagocyte_y + dy_phagocyte, 1), grid_size)
    }

    # ========================================================================
    # STEP 4: ABM-PDE INTERACTIONS
    # ========================================================================

    # 4.1 Decay macrophage memory (exponential decay)
    phagocyte_pathogen_memory = phagocyte_pathogen_memory * memory_decay_rate
    phagocyte_commensal_memory = phagocyte_commensal_memory * memory_decay_rate

    # 4.2 Macrophage engulfment (samples from PDE, removes from PDE)
    for (i in 1:n_phagocytes) {
      px = phagocyte_x[i]
      py = phagocyte_y[i]

      # Sample pathogen engulfment from local concentration
      local_P = P_field[py, px]
      engulf_prob_P = phagocyte_activity_engulf[i]
      n_pathogens_engulfed = rpois(1, engulf_prob_P * local_P)
      n_pathogens_engulfed = min(n_pathogens_engulfed, local_P)  # Can't engulf more than present

      if (n_pathogens_engulfed > 0) {
        P_field[py, px] = P_field[py, px] - n_pathogens_engulfed
        phagocyte_pathogen_memory[i] = phagocyte_pathogen_memory[i] + n_pathogens_engulfed
        pathogens_killed_by_Mac[phagocyte_phenotype[i] + 1] =
          pathogens_killed_by_Mac[phagocyte_phenotype[i] + 1] + n_pathogens_engulfed
      }

      # Sample commensal engulfment from local concentration
      local_C = C_field[py, px]
      engulf_prob_C = phagocyte_activity_engulf[i]
      n_commensals_engulfed = rpois(1, engulf_prob_C * local_C)
      n_commensals_engulfed = min(n_commensals_engulfed, local_C)

      if (n_commensals_engulfed > 0) {
        C_field[py, px] = C_field[py, px] - n_commensals_engulfed
        phagocyte_commensal_memory[i] = phagocyte_commensal_memory[i] + n_commensals_engulfed
        commensals_killed_by_Mac[phagocyte_phenotype[i] + 1] =
          commensals_killed_by_Mac[phagocyte_phenotype[i] + 1] + n_commensals_engulfed
      }
    }

    # 4.3 ROS killing (PDE × PDE interaction)
    # Kill rate proportional to ROS concentration × microbe concentration
    P_killed_by_ROS = pmin(P_field, th_ROS_microbe * ROS * P_field)
    C_killed_by_ROS = pmin(C_field, th_ROS_microbe * ROS * C_field)

    P_field = P_field - P_killed_by_ROS
    C_field = C_field - C_killed_by_ROS

    pathogens_killed_by_ROS = pathogens_killed_by_ROS + sum(P_killed_by_ROS)
    commensals_killed_by_ROS = commensals_killed_by_ROS + sum(C_killed_by_ROS)

    # ========================================================================
    # STEP 5: MACROPHAGE POLARIZATION (uses decaying memory!)
    # ========================================================================

    M0_indices = which(phagocyte_phenotype == 0)
    M1_indices = which(phagocyte_phenotype == 1)
    M2_indices = which(phagocyte_phenotype == 2)

    # Process M0 phagocytes (candidates for activation)
    if (length(M0_indices) > 0) {
      signals = calculate_phagocyte_signals_cpp(
        M0_indices, phagocyte_x, phagocyte_y,
        act_radius_DAMPs, act_radius_SAMPs, DAMPs, SAMPs, grid_size
      )
      avg_DAMPs_vec = signals$avg_DAMPs
      avg_SAMPs_vec = signals$avg_SAMPs

      for (idx in seq_along(M0_indices)) {
        i = M0_indices[idx]
        avg_DAMPs = avg_DAMPs_vec[idx]
        avg_SAMPs = avg_SAMPs_vec[idx]

        if (avg_DAMPs >= activation_threshold_DAMPs && avg_DAMPs > avg_SAMPs) {
          phagocyte_phenotype[i] = 1
          phagocyte_active_age[i] = 1
          phagocyte_activity_ROS[i] = activity_ROS_M1
          phagocyte_activity_engulf[i] = activity_engulf_M1
        } else if (avg_SAMPs >= activation_threshold_SAMPs && avg_SAMPs > avg_DAMPs) {
          phagocyte_phenotype[i] = 2
          phagocyte_active_age[i] = 1
          phagocyte_activity_ROS[i] = activity_ROS_M2
          phagocyte_activity_engulf[i] = activity_engulf_M2
        }
      }
    }

    # Process M1/M2 phagocytes (can switch or deactivate)
    active_indices = c(M1_indices, M2_indices)
    if (length(active_indices) > 0) {
      phagocyte_active_age[active_indices] = phagocyte_active_age[active_indices] + 1
      old_enough = phagocyte_active_age[active_indices] >= active_age_limit
      candidates = active_indices[old_enough]

      if (length(candidates) > 0) {
        signals = calculate_phagocyte_signals_cpp(
          candidates, phagocyte_x, phagocyte_y,
          act_radius_DAMPs, act_radius_SAMPs, DAMPs, SAMPs, grid_size
        )
        avg_DAMPs_vec = signals$avg_DAMPs
        avg_SAMPs_vec = signals$avg_SAMPs

        for (idx in seq_along(candidates)) {
          i = candidates[idx]
          avg_DAMPs = avg_DAMPs_vec[idx]
          avg_SAMPs = avg_SAMPs_vec[idx]

          if(macspec_on==1){
            # Macrophage specificity: use MEMORY to inform polarization
            total_memory = phagocyte_pathogen_memory[i] + phagocyte_commensal_memory[i]

            pathogen_memory_dominant = FALSE
            commensal_memory_dominant = FALSE

            if (total_memory > 0) {
              rat_com_pat_real = phagocyte_commensal_memory[i] / total_memory
              alpha = (1 - mac_discrimination_efficiency) * 1 +
                mac_discrimination_efficiency * (rat_com_pat_real * precision_mac)
              beta = (1 - mac_discrimination_efficiency) * 1 +
                mac_discrimination_efficiency * ((1 - rat_com_pat_real) * precision_mac)

              rat_com_pat = sample_rbeta(alpha, beta)

              pathogen_memory_dominant  = rat_com_pat <= (1 - mac_rat_com_pat_threshold)
              commensal_memory_dominant = (rat_com_pat > mac_rat_com_pat_threshold)
            }

            DAMPs_dominant = (avg_DAMPs >= activation_threshold_DAMPs && avg_DAMPs > avg_SAMPs)
            SAMPs_dominant = (avg_SAMPs >= activation_threshold_SAMPs && avg_SAMPs > avg_DAMPs)

            # POLARIZATION LOGIC: Danger dominates
            if (DAMPs_dominant || pathogen_memory_dominant) {
              phagocyte_phenotype[i] = 1
              phagocyte_active_age[i] = 1
              phagocyte_activity_ROS[i] = activity_ROS_M1
              phagocyte_activity_engulf[i] = activity_engulf_M1
            } else if (SAMPs_dominant && commensal_memory_dominant) {
              phagocyte_phenotype[i] = 2
              phagocyte_active_age[i] = 1
              phagocyte_activity_ROS[i] = activity_ROS_M2
              phagocyte_activity_engulf[i] = activity_engulf_M2
            } else if (avg_SAMPs < activation_threshold_SAMPs && avg_DAMPs < activation_threshold_DAMPs) {
              phagocyte_phenotype[i] = 0
              phagocyte_active_age[i] = 0
              phagocyte_activity_ROS[i] = activity_ROS_M0
              phagocyte_activity_engulf[i] = activity_engulf_M0
            }
          }else{
            # Vanilla mode: only environmental signals
            if (avg_DAMPs >= activation_threshold_DAMPs && avg_DAMPs > avg_SAMPs) {
              phagocyte_phenotype[i] = 1
              phagocyte_active_age[i] = 1
              phagocyte_activity_ROS[i] = activity_ROS_M1
              phagocyte_activity_engulf[i] = activity_engulf_M1
            } else if (avg_SAMPs >= activation_threshold_SAMPs && avg_SAMPs > avg_DAMPs) {
              phagocyte_phenotype[i] = 2
              phagocyte_active_age[i] = 1
              phagocyte_activity_ROS[i] = activity_ROS_M2
              phagocyte_activity_engulf[i] = activity_engulf_M2
            } else if (avg_SAMPs < activation_threshold_SAMPs && avg_DAMPs < activation_threshold_DAMPs) {
              phagocyte_phenotype[i] = 0
              phagocyte_active_age[i] = 0
              phagocyte_activity_ROS[i] = activity_ROS_M0
              phagocyte_activity_engulf[i] = activity_engulf_M0
            }
          }
        }
      }
    }

    # ========================================================================
    # STEP 6: TREG ACTIVATION (by nearby macrophages)
    # ========================================================================

    M1_phagocyte_indices = which(phagocyte_phenotype == 1)
    M2_phagocyte_indices = which(phagocyte_phenotype == 2)
    M_activate_phagocyte_indices = c(M1_phagocyte_indices, M2_phagocyte_indices)

    if (length(M_activate_phagocyte_indices) > 0) {
      for (i in M_activate_phagocyte_indices) {
        px = phagocyte_x[i]
        py = phagocyte_y[i]

        nearby_treg_indices = find_nearby_tregs_cpp(
          px, py, treg_x, treg_y, treg_vicinity_effect
        )

        if (length(nearby_treg_indices) > 0) {
          total_memory = phagocyte_pathogen_memory[i] + phagocyte_commensal_memory[i]

          if (total_memory > 0) {
            rat_com_pat_real = phagocyte_commensal_memory[i] / total_memory
            alpha = (1 - treg_discrimination_efficiency) * 1 +
              treg_discrimination_efficiency * (rat_com_pat_real * precision_treg)
            beta = (1 - treg_discrimination_efficiency) * 1 +
              treg_discrimination_efficiency * ((1 - rat_com_pat_real) * precision_treg)

            rat_com_pat = sample_rbeta(alpha, beta)

            if (rat_com_pat > rat_com_pat_threshold) {
              treg_phenotype[nearby_treg_indices] = 1
              treg_activity_SAMPs_binary[nearby_treg_indices] = 1
              treg_active_age[nearby_treg_indices] = 1
            }
          }
        }
      }
    }

    # ========================================================================
    # STEP 7: TREG AGING & DEACTIVATION
    # ========================================================================

    active_treg_indices = which(treg_phenotype == 1)
    if (length(active_treg_indices) > 0) {
      old_tregs = active_treg_indices[treg_active_age[active_treg_indices] >= active_age_limit]
      young_tregs = active_treg_indices[treg_active_age[active_treg_indices] < active_age_limit]

      if (length(young_tregs) > 0) {
        treg_active_age[young_tregs] = treg_active_age[young_tregs] + 1
      }

      if (length(old_tregs) > 0) {
        treg_phenotype[old_tregs] = 0
        treg_active_age[old_tregs] = 0
        treg_activity_SAMPs_binary[old_tregs] = 0
      }
    }

    # ========================================================================
    # STEP 8: EPITHELIAL DYNAMICS
    # ========================================================================

    # 8.1 Calculate ROS at epithelium
    epithelium_x = epithelium$x
    ros_means = calculate_epithelial_ros_cpp(
      epithelium_x, act_radius_ROS, ROS, grid_size
    )

    # 8.2 Injury from pathogens at epithelium
    for (k in 1:nrow(epithelium)) {
      x = epithelium$x[k]
      pathogen_load = P_field[1, x]
      epithelium$level_injury[k] = epithelium$level_injury[k] +
        logistic_scaled_0_to_5_quantized(pathogen_load)
    }

    # 8.3 Injury from ROS
    epithelium$level_injury = epithelium$level_injury + as.integer(ros_means > th_ROS_epith_recover)

    # 8.4 Apply maximum injury constraint
    epithelium$level_injury = pmin(epithelium$level_injury, max_level_injury)

    # 8.5 Stochastic recovery
    for (k in 1:nrow(epithelium)) {
      if (epithelium$level_injury[k] > 0 && runif(1) < epith_recovery_chance) {
        epithelium$level_injury[k] = max(0, epithelium$level_injury[k] - 1)
      }
    }

    # ========================================================================
    # STEP 9: RECORD LONGITUDINAL DATA
    # ========================================================================

    epithelium_longitudinal[t, ] = as.numeric(table(factor(epithelium$level_injury, levels = 0:5)))

    phagocyte_counts = c(
      sum(phagocyte_phenotype == 0),
      sum(phagocyte_phenotype == 1),
      sum(phagocyte_phenotype == 2)
    )
    macrophages_longitudinal[t, ] = phagocyte_counts

    # Total microbe counts (sum of concentration fields)
    microbes_longitudinal[t, ] = c(sum(C_field), sum(P_field))

    tregs_longitudinal[t, ] = c(sum(treg_phenotype == 0), sum(treg_phenotype == 1))

    microbes_cumdeath_longitudinal[t, ] = c(
      commensals_killed_by_ROS, commensals_killed_by_Mac,
      pathogens_killed_by_ROS, pathogens_killed_by_Mac
    )
  }

  # ==========================================================================
  # SAVE LONGITUDINAL DATA
  # ==========================================================================

  longitudinal_df = data.frame(
    epithelium_longitudinal,
    macrophages_longitudinal,
    microbes_longitudinal,
    tregs_longitudinal,
    microbes_cumdeath_longitudinal
  )

  colnames(longitudinal_df) = colnames_insert_1level

  longitudinal_df$t = 1:t_max
  longitudinal_df$control = control
  longitudinal_df$sterile = sterile
  longitudinal_df$macspec_on = macspec_on
  longitudinal_df$tregs_on = allow_tregs
  longitudinal_df$randomize_tregs = randomize_tregs
  longitudinal_df$param_set_id = param_set_use$param_set_id
  longitudinal_df$rep_id = reps_in

  longitudinal_df = longitudinal_df %>% dplyr::mutate(epithelial_score = 6*epithelial_healthy+
                                                        5*epithelial_inj_1+
                                                        4*epithelial_inj_2+
                                                        3*epithelial_inj_3+
                                                        2*epithelial_inj_4+
                                                        1*epithelial_inj_5)

  longitudinal_df$time_ss = steady_state_idx(longitudinal_df$epithelial_score)

  longitudinal_df = longitudinal_df %>%
    dplyr::select(t,
                  control,
                  sterile,
                  tregs_on,
                  macspec_on,
                  randomize_tregs,
                  param_set_id,
                  rep_id,
                  epithelial_score,
                  time_ss,
                  everything())

  longitudinal_df_keep = rbind(longitudinal_df_keep, longitudinal_df)
}
