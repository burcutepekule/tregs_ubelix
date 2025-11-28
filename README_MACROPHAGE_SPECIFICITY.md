# Macrophage Specificity Implementation

## Overview

This implementation adds **macrophage discrimination** capability to the simulation, allowing macrophages to use both environmental signals (DAMPs/SAMPs) AND engulfment patterns (pathogens vs commensals) to determine their polarization state.

## Biological Rationale

**Danger signals dominate immune responses:**
- If there are high DAMPs OR high pathogen engulfment → M1 (pro-inflammatory)
- M2 (anti-inflammatory) polarization requires BOTH:
  - High SAMPs in environment
  - Predominantly commensal engulfment
- This reflects the immune system's bias toward safety: when in doubt, respond to danger

## New Parameters

### 1. `mac_discrimination_efficiency` (range: 0-1)
- How well macrophages can discriminate between pathogens and commensals they engulf
- **0 = random** (cannot distinguish, sees 50/50 regardless of actual ratio)
- **1 = perfect** (perfectly perceives actual pathogen/commensal ratio)
- Uses same beta distribution approach as Treg discrimination

### 2. `mac_rat_com_pat_threshold` (range: 0.5-1)
- Threshold commensal ratio needed for M2 polarization
- **0.5** = requires equal or more commensals than pathogens
- **1.0** = requires 100% commensals (no pathogens)
- Higher values make M2 polarization more stringent

## Polarization Logic

### M1 Polarization (Pro-inflammatory)
Occurs if **EITHER** condition is true:
1. **Environmental danger**: `(DAMPs > threshold) AND (DAMPs > SAMPs)`
2. **Pathogen engulfment dominance**: `pathogen_ratio > (1 - mac_rat_com_pat_threshold)`

### M2 Polarization (Anti-inflammatory)
Occurs if **BOTH** conditions are true:
1. **Environmental safety**: `(SAMPs > threshold) AND (SAMPs > DAMPs)`
2. **Commensal engulfment dominance**: `commensal_ratio > mac_rat_com_pat_threshold`

### M0 (Resting)
- If neither M1 nor M2 conditions are met
- OR if both environmental signals are below threshold

## Implementation Details

### Engulfment Ratio Calculation
Similar to Treg activation, uses beta distribution to model discrimination noise:

```r
rat_com_pat_real = num_commensals / (num_commensals + num_pathogens)

alpha = (1 - mac_discrimination_efficiency) * 1 +
        mac_discrimination_efficiency * (rat_com_pat_real * precision_mac)

beta = (1 - mac_discrimination_efficiency) * 1 +
       mac_discrimination_efficiency * ((1 - rat_com_pat_real) * precision_mac)

rat_com_pat = sample_rbeta(alpha, beta)  # Perceived ratio (with noise)
```

Where: `precision_mac = 10 * exp(5 * mac_discrimination_efficiency)`

### Key Decision Points

The combinatorial logic handles four possible scenarios:

| Environment | Engulfment | Result | Rationale |
|-------------|------------|--------|-----------|
| High DAMPs  | High pathogen | M1 | ✓ Concordant danger |
| High DAMPs  | High commensal | M1 | ✓ Environmental danger dominates |
| High SAMPs  | High pathogen | M1 | ✓ Engulfed danger dominates |
| High SAMPs  | High commensal | M2 | ✓ Concordant safety |

## Files Created

### 1. **UBX_paramgen_lhs_macspec.R**
- Generates LHS parameter samples including new macrophage parameters
- Output: `lhs_parameters_macspec_ubelix.csv`

### 2. **MISC/ASSIGN_PARAMETERS_MACSPEC.R**
- Loads parameters from CSV including:
  - `mac_discrimination_efficiency`
  - `mac_rat_com_pat_threshold`
- Calculates `precision_mac` for beta distribution sampling

### 3. **MISC/RUN_REPS_CPP_MACSPEC.R**
- Modified simulation loop with new polarization logic
- Changes in two sections:
  - **M0 phagocyte processing** (lines ~380-429): Initial activation
  - **M1/M2 phagocyte processing** (lines ~462-515): Repolarization/reversion

### 4. **UBX_datagen_cpp_macspec.R**
- Main data generation script
- Uses macrophage specificity parameter file
- Output directory: `mass_sim_results_R_cpp_macspec/`

## Usage

### 1. Generate Parameters
```bash
Rscript UBX_paramgen_lhs_macspec.R
```
Creates: `lhs_parameters_macspec_ubelix.csv` with 100,000 parameter sets

### 2. Run Simulations
```bash
# Run chunk n2 out of n1 total chunks
Rscript UBX_datagen_cpp_macspec.R <n1> <n2>

# Example: Run chunk 1 out of 100
Rscript UBX_datagen_cpp_macspec.R 100 1
```

### 3. Output
Results saved to: `/storage/homefs/bt25p365/tregs/mass_sim_results_R_cpp_macspec/`

Files named: `longitudinal_df_param_set_id_{id}_sterile_{0|1}_tregs_{0|1}_trnd_{0|1}.rds`

## Comparison with Original Implementation

| Aspect | Original | With Macrophage Specificity |
|--------|----------|----------------------------|
| Polarization cues | DAMPs vs SAMPs only | DAMPs/SAMPs + Engulfment pattern |
| M1 trigger | DAMPs > SAMPs | DAMPs > SAMPs OR pathogen engulfment |
| M2 trigger | SAMPs > DAMPs | SAMPs > DAMPs AND commensal engulfment |
| Parameters | 23 | 25 (+2 mac specificity) |
| Discrimination | Tregs only | Tregs + Macrophages |

## Expected Biological Impacts

### Higher `mac_discrimination_efficiency`:
- More accurate M1/M2 decisions based on actual microbe composition
- Stronger differentiation between pathogenic and commensal environments
- May improve outcomes in mixed microbiome scenarios

### Higher `mac_rat_com_pat_threshold`:
- More stringent M2 polarization requirements
- Could reduce inappropriate anti-inflammatory responses
- May increase M1 bias overall

### Low discrimination + stringent threshold:
- Difficult to achieve M2 (noisy perception + high requirement)
- Likely to increase M1 polarization
- Could lead to excessive inflammation

## Random Number Stream Synchronization

**IMPORTANT**: The implementation maintains random number stream synchronization by:
1. ALWAYS sampling from beta distribution when macrophages have engulfed bacteria
2. This occurs even if the result doesn't affect polarization
3. Ensures reproducibility across parameter variations

## Testing Recommendations

1. **Sanity checks:**
   - Pure pathogen environment → mostly M1?
   - Pure commensal environment → mostly M2?
   - Mixed environments → intermediate?

2. **Parameter sensitivity:**
   - Does `mac_discrimination_efficiency = 0` behave randomly?
   - Does `mac_discrimination_efficiency = 1` show clear separation?
   - Effect of `mac_rat_com_pat_threshold` on M2 abundance?

3. **Comparison studies:**
   - Original vs macspec with same base parameters
   - How much does engulfment pattern matter relative to environmental signals?

## Notes

- Both macrophages and Tregs now have independent discrimination parameters
- Maintains C++ acceleration compatibility
- Backward compatible: can generate data without macrophage specificity using original scripts
- Consider: Do macrophages and Tregs see the same "commensal" definition?
