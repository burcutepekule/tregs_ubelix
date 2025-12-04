# Hybrid PDE-ABM Model: Complete Algorithm

## Architecture

### PDE Components (continuous fields on 25×25 grid)
- `P_field[i,j]` - pathogen concentration
- `C_field[i,j]` - commensal concentration
- `DAMPs[i,j]` - DAMP concentration
- `SAMPs[i,j]` - SAMP concentration
- `ROS[i,j]` - ROS concentration

### ABM Components (discrete agents)

**Epithelial cells** (n=25, fixed at y=1):
- `epithelium_x[k]` - x position (1 to 25)
- `epithelium_injury[k]` - injury level (0 to 5)

**Macrophages** (n=~31):
- `mac_x[i], mac_y[i]` - position
- `mac_phenotype[i]` - {0=M0, 1=M1, 2=M2}
- `mac_pathogen_memory[i]` - decaying memory of pathogen engulfment
- `mac_commensal_memory[i]` - decaying memory of commensal engulfment
- `mac_active_age[i]`
- `mac_activity_ROS[i]`
- `mac_activity_engulf[i]`

**Tregs** (n=~31):
- `treg_x[i], treg_y[i]` - position
- `treg_phenotype[i]` - {0=resting, 1=active}
- `treg_active_age[i]`
- `treg_activity_SAMPs[i]`

---

## Timestep Algorithm

### **STEP 1: Agent Actions → PDE Sources**

#### 1.1 SAMPs production (from active Tregs)
```r
for (i in active_tregs) {
  x = treg_x[i]
  y = treg_y[i]
  S[y, x] = S[y, x] + treg_activity_SAMPs[i] * add_SAMPs * allow_tregs
}
```

#### 1.2 ROS production (from M1 macrophages)
```r
for (i in M1_macs) {
  x = mac_x[i]
  y = mac_y[i]
  R[y, x] = R[y, x] + mac_activity_ROS[i] * add_ROS  # CONTROL: add_ROS=0
}
```

#### 1.3 DAMP production (from epithelial injury + microbe contact)
```r
# From epithelial injury
for (k in 1:n_epithelium) {
  x = epithelium_x[k]
  D[1, x] = D[1, x] + epithelium_injury[k] * add_DAMPs
}

# From microbe density at epithelium
for (k in 1:n_epithelium) {
  x = epithelium_x[k]
  microbe_density = P[1, x] + C[1, x]
  D[1, x] = D[1, x] + logistic_scaled(microbe_density) * add_DAMPs
}
```

---

### **STEP 2: PDE Evolution**

#### 2.1 Pathogen/Commensal entry (from epithelium)
```r
# Pathogens leak through injured epithelium
for (k in 1:n_epithelium) {
  x = epithelium_x[k]
  inj = epithelium_injury[k]
  P[1, x] = P[1, x] + rate_leak_pathogen * inj * dt
}

# Commensals leak through baseline + injury
for (k in 1:n_epithelium) {
  x = epithelium_x[k]
  inj = epithelium_injury[k]
  C[1, x] = C[1, x] + (rate_leak_commensal_baseline + rate_leak_commensal_injury * inj) * dt
}
```

#### 2.2 Diffusion step (all fields)
```r
P = diffuse_matrix(P, D_pathogen, max_value)
C = diffuse_matrix(C, D_commensal, max_value)
D = diffuse_matrix(D, D_damp, max_value)
S = diffuse_matrix(S, D_samp, max_value)
R = diffuse_matrix(R, D_ros, max_value)
```

#### 2.3 Decay
```r
D = D - decay_damp * D
S = S - decay_samp * S
R = R - decay_ros * R
# Note: P, C don't decay naturally (they're bacteria, not chemicals)
```

---

### **STEP 3: ABM Movement**

#### 3.1 Macrophage chemotaxis (toward DAMPs)
```r
for (i in 1:n_macrophages) {
  x = mac_x[i]
  y = mac_y[i]

  # Calculate DAMP gradient in 8-neighborhood
  neighbors = get_8_neighbors(x, y)
  damp_values = D[neighbors]

  # Move probabilistically toward highest DAMP
  probs = damp_values / sum(damp_values)
  chosen = sample(neighbors, prob=probs)
  mac_x[i] = chosen$x
  mac_y[i] = chosen$y
}
```

#### 3.2 Treg chemotaxis (toward DAMPs, unless randomized)
```r
if (randomize_tregs == 0) {
  # Same as macrophages
  for (i in 1:n_tregs) {
    # ... move toward DAMP gradient
  }
} else {
  # Random walk
  for (i in 1:n_tregs) {
    treg_x[i] = treg_x[i] + sample(c(-1,0,1), 1)
    treg_y[i] = treg_y[i] + sample(c(-1,0,1), 1)
    # Clamp to boundaries
  }
}
```

---

### **STEP 4: ABM Interactions (Agents ↔ PDE fields)**

#### 4.1 Macrophage engulfment (removes from PDE fields)
```r
for (i in 1:n_macrophages) {
  x = mac_x[i]
  y = mac_y[i]

  # Sample local pathogen concentration
  engulf_rate_P = mac_activity_engulf[i] * P[y, x] * dt
  n_pathogens_engulfed = rpois(1, lambda = engulf_rate_P)

  # Sample local commensal concentration
  engulf_rate_C = mac_activity_engulf[i] * C[y, x] * dt
  n_commensals_engulfed = rpois(1, lambda = engulf_rate_C)

  # Update PDE fields (remove engulfed microbes)
  P[y, x] = max(0, P[y, x] - n_pathogens_engulfed)
  C[y, x] = max(0, C[y, x] - n_commensals_engulfed)

  # Update macrophage memory (with new engulfment)
  mac_pathogen_memory[i] = mac_pathogen_memory[i] + n_pathogens_engulfed
  mac_commensal_memory[i] = mac_commensal_memory[i] + n_commensals_engulfed
}
```

#### 4.2 ROS killing (PDE removes from PDE)
```r
# Pathogens killed where ROS is high
kill_rate_P = th_ROS_microbe * R * P * dt
P = P - kill_rate_P
P = pmax(P, 0)  # ensure non-negative

# Commensals killed where ROS is high
kill_rate_C = th_ROS_microbe * R * C * dt
C = C - kill_rate_C
C = pmax(C, 0)
```

---

### **STEP 5: Macrophage Polarization (uses decaying memory!)**

#### 5.1 Decay existing memory
```r
memory_half_life = 20  # timesteps (PARAMETER!)
decay_rate = exp(-log(2) / memory_half_life)

mac_pathogen_memory = mac_pathogen_memory * decay_rate
mac_commensal_memory = mac_commensal_memory * decay_rate
```

#### 5.2 Sense local signals
```r
for (i in 1:n_macrophages) {
  x = mac_x[i]
  y = mac_y[i]

  # Average DAMP/SAMP in neighborhood
  avg_damp = mean(D[y + (-radius:radius), x + (-radius:radius)])
  avg_samp = mean(S[y + (-radius:radius), x + (-radius:radius)])

  # Store for polarization decision
  mac_local_damp[i] = avg_damp
  mac_local_samp[i] = avg_samp
}
```

#### 5.3 Polarization logic

**Vanilla mode (macspec_on=0):**
```r
for (i in M0_macrophages) {
  if (avg_damp > threshold_damp && avg_damp > avg_samp) {
    mac_phenotype[i] = 1  # M1
    mac_activity_ROS[i] = activity_ROS_M1_baseline
    mac_activity_engulf[i] = activity_engulf_M1_baseline
  } else if (avg_samp > threshold_samp && avg_samp > avg_damp) {
    mac_phenotype[i] = 2  # M2
    mac_activity_ROS[i] = activity_ROS_M2_baseline
    mac_activity_engulf[i] = activity_engulf_M2_baseline
  }
}
```

**Macrophage specificity mode (macspec_on=1):**
```r
for (i in candidate_macrophages) {
  # Calculate commensal ratio from MEMORY
  total_memory = mac_pathogen_memory[i] + mac_commensal_memory[i]

  if (total_memory > 0) {
    # Sample with discrimination efficiency
    rat_com_pat_real = mac_commensal_memory[i] / total_memory

    alpha = (1 - mac_discrimination_efficiency) * 1 +
            mac_discrimination_efficiency * (rat_com_pat_real * precision_mac)
    beta = (1 - mac_discrimination_efficiency) * 1 +
           mac_discrimination_efficiency * ((1 - rat_com_pat_real) * precision_mac)

    rat_com_pat = rbeta(1, alpha, beta)

    pathogen_memory_dominant = (rat_com_pat <= (1 - mac_rat_com_pat_threshold))
    commensal_memory_dominant = (rat_com_pat > mac_rat_com_pat_threshold)
  } else {
    pathogen_memory_dominant = FALSE
    commensal_memory_dominant = FALSE
  }

  # Determine environmental signals
  damp_dominant = (avg_damp >= threshold_damp && avg_damp > avg_samp)
  samp_dominant = (avg_samp >= threshold_samp && avg_samp > avg_damp)

  # POLARIZATION LOGIC: Danger dominates
  if (damp_dominant || pathogen_memory_dominant) {
    # M1: Either environmental danger OR pathogen memory
    mac_phenotype[i] = 1
    mac_active_age[i] = 1
    mac_activity_ROS[i] = activity_ROS_M1_baseline
    mac_activity_engulf[i] = activity_engulf_M1_baseline

  } else if (samp_dominant && commensal_memory_dominant) {
    # M2: Both environmental safety AND commensal memory required
    mac_phenotype[i] = 2
    mac_active_age[i] = 1
    mac_activity_ROS[i] = activity_ROS_M2_baseline
    mac_activity_engulf[i] = activity_engulf_M2_baseline

  } else if (avg_samp < threshold_samp && avg_damp < threshold_damp) {
    # Revert to M0 if both signals low
    mac_phenotype[i] = 0
    mac_active_age[i] = 0
    mac_activity_ROS[i] = activity_ROS_M0_baseline
    mac_activity_engulf[i] = activity_engulf_M0_baseline
  }
}
```

---

### **STEP 6: Treg Activation**

```r
for (i in M1_or_M2_macrophages) {
  # Find nearby Tregs
  nearby_tregs = find_tregs_within_radius(mac_x[i], mac_y[i], treg_vicinity_effect)

  if (length(nearby_tregs) > 0) {
    # Sample antigen presentation using macrophage's MEMORY
    total_memory = mac_pathogen_memory[i] + mac_commensal_memory[i]

    if (total_memory > 0) {
      rat_com_pat_real = mac_commensal_memory[i] / total_memory

      alpha = (1 - treg_discrimination_efficiency) * 1 +
              treg_discrimination_efficiency * (rat_com_pat_real * precision_treg)
      beta = (1 - treg_discrimination_efficiency) * 1 +
             treg_discrimination_efficiency * ((1 - rat_com_pat_real) * precision_treg)

      rat_com_pat = rbeta(1, alpha, beta)

      # If commensal-biased, activate nearby Tregs
      if (rat_com_pat > rat_com_pat_threshold) {
        treg_phenotype[nearby_tregs] = 1
        treg_activity_SAMPs[nearby_tregs] = 1
        treg_active_age[nearby_tregs] = 1
      }
    }
  }
}
```

---

### **STEP 7: Epithelial Dynamics**

#### 7.1 Injury from pathogens
```r
for (k in 1:n_epithelium) {
  x = epithelium_x[k]
  pathogen_load = P[1, x]

  injury_increase = logistic_scaled(pathogen_load)
  epithelium_injury[k] = epithelium_injury[k] + injury_increase
}
```

#### 7.2 Injury from ROS
```r
for (k in 1:n_epithelium) {
  x = epithelium_x[k]
  local_ros = mean(R[1, max(1, x-radius):min(grid_size, x+radius)])

  if (local_ros > th_ROS_epith_recover) {
    epithelium_injury[k] = epithelium_injury[k] + 1
  }
}
```

#### 7.3 Recovery
```r
for (k in 1:n_epithelium) {
  if (epithelium_injury[k] > 0 && runif(1) < epith_recovery_chance) {
    epithelium_injury[k] = epithelium_injury[k] - 1
  }
}

# Clamp to max injury
epithelium_injury = pmin(epithelium_injury, max_level_injury)
```

---

### **STEP 8: Aging & Deactivation**

```r
# Age active macrophages
mac_active_age[active_macs] = mac_active_age[active_macs] + 1

# Deactivate old macrophages (handle in STEP 5)

# Age active Tregs
treg_active_age[active_tregs] = treg_active_age[active_tregs] + 1

# Deactivate old Tregs
old_tregs = which(treg_phenotype == 1 & treg_active_age >= active_age_limit)
treg_phenotype[old_tregs] = 0
treg_active_age[old_tregs] = 0
treg_activity_SAMPs[old_tregs] = 0
```

---

## Key Parameters to Add

```r
# Memory decay
memory_half_life = 20  # How fast engulfment memory decays

# PDE diffusion rates (need to calibrate!)
D_pathogen = 0.1    # Pathogens diffuse like particles
D_commensal = 0.1   # Commensals diffuse like particles
```

---

## Expected Speedup

**Old ABM:**
- Track every individual microbe (could be 1000s)
- Random walks for each
- Collision detection for each

**New Hybrid:**
- Microbes as continuous field (fast PDE solve)
- Only ~60 agents to track
- **Expected: 50-200× speedup**

---

## Implementation Notes

1. **Microbe units**: PDE fields use concentrations (arbitrary units). Calibrate so `P[i,j]~1.0` ≈ 1 bacterium in old model
2. **Engulfment stochasticity**: Use Poisson sampling from concentration field
3. **Memory tuning**: `memory_half_life` is critical parameter - test 10-50 timesteps
4. **Spatial resolution**: Can reduce grid to 15×15 for even more speed

---

## Control Scenario Implementation

Same as before - just set:
```r
if (control == 1) {
  add_ROS = 0
  activity_ROS_M1_baseline = 0
}
```

This disables ROS production but keeps engulfment functional.
