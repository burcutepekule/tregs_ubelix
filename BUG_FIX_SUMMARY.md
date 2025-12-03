# Bug Fix: C++ vs R Implementation Discrepancy

## Problem
When running `UBX_datagen_vanilla_local_cpp.R` with `cpp_on = T` vs `cpp_on = F` (and `one_level = F`), different results were obtained despite using the same random seed.

## Root Cause
The bug was in the **epithelial ROS calculation** where the R fallback implementation created duplicate indices at grid boundaries, while the C++ implementation did not.

### R Version (Buggy):
```r
x_coordinates = pmax(1, pmin(grid_size, (px - act_radius_ROS):(px + act_radius_ROS)))
ros_means[i] = mean(ROS[1, x_coordinates])
```

**Problem**: `pmax` and `pmin` operate element-wise on the range, which can create duplicates.

**Example** (px=1, act_radius_ROS=1, grid_size=25):
1. Range: `(1-1):(1+1)` = `0:2` = `c(0, 1, 2)`
2. `pmin(25, c(0,1,2))` = `c(0, 1, 2)`
3. `pmax(1, c(0,1,2))` = `c(1, 1, 2)` ← **DUPLICATE!**
4. `mean(ROS[1, c(1,1,2)])` counts `ROS[1,1]` **TWICE**

### C++ Version (Correct):
```cpp
int x_start = std::max(0, x0 - act_radius_ROS);
int x_end = std::min(grid_size - 1, x0 + act_radius_ROS);
for (int xi = x_start; xi <= x_end; xi++) {
  sum += ROS(0, xi);
  count++;
}
```

**Correct behavior**: Clamps the range endpoints, never creates duplicates.

## Impact
- Epithelial cells at grid boundaries (x=1 and x=grid_size) received incorrect ROS exposure calculations in the R version
- This affected epithelial injury progression differently between C++ and R implementations
- The difference propagated through the simulation, causing divergent outcomes

## Fix Applied
Updated the R fallback code in all simulation files to match the C++ implementation:

```r
# FIXED: Clamp range endpoints to avoid duplicate indices at boundaries
x_start = max(1, px - act_radius_ROS)
x_end = min(grid_size, px + act_radius_ROS)
x_coordinates = x_start:x_end
ros_means[i] = mean(ROS[1, x_coordinates])
```

### Files Modified:
1. `MISC/RUN_REPS_CPP.R` (line 666-675)
2. `MISC/RUN_REPS.R` (line 540-545)
3. `MISC/RUN_REPS_OPTIMIZED.R` (line 618-623)
4. `MISC/RUN_REPS_CPP_ONELEVEL.R` (line 666-674)
5. `MISC/RUN_REPS_CPP_MACSPEC_AND_VANILLA_ONELEVEL.R` (line 680-688)
6. `MISC/RUN_REPS_CPP_MACSPEC_AND_VANILLA.R` (line 715-723)
7. `MISC/RUN_REPS_CPP_MACSPEC.R` (line 728-736)
8. `x_UBX_datagen_v0.R` (line 700-705)

## Verification
After this fix:
- C++ implementation (cpp_on = T) and R implementation (cpp_on = F) should now produce identical results
- No more duplicate index calculations at boundaries
- Consistent epithelial ROS exposure calculations across both implementations

## Note
The C++ implementation was correct all along. The bug was only in the R fallback version that gets used when `cpp_on = F` or when C++ compilation fails.
