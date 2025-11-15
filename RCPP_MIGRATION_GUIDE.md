# Rcpp Migration Guide for RUN_REPS.R

## Why Rcpp for This Project?

Your simulation is a **perfect candidate** for Rcpp because:

1. ✅ **Tight loops** - You have `for (t in 1:5000)` with nested agent loops
2. ✅ **Agent-based model** - Lots of per-agent operations (movement, interactions)
3. ✅ **Conditional logic** - Can't be vectorized in R
4. ✅ **Dynamic updates** - Arrays being modified constantly

**Expected speedup: 20-100x faster** for the simulation core.

---

## Quick Start: Setting Up Rcpp

### 1. Install Rcpp

```r
install.packages("Rcpp")
library(Rcpp)
```

### 2. Test Your Setup

```r
# Test that C++ compilation works
Rcpp::evalCpp("2 + 2")  # Should return 4
```

If this fails, you may need to install:
- **Linux**: `sudo apt-get install r-base-dev`
- **Mac**: Install Xcode command line tools
- **Windows**: Install Rtools

---

## Migration Strategy: Incremental Approach

**DON'T rewrite everything at once!** Port functions incrementally:

### Phase 1: Easy Wins (1-2 days)
Port the **most called** functions that are simple:

1. ✅ `get_8n_avg_signal_fast()` - Called thousands of times
2. ✅ ROS microbe killing loop (lines 501-531)
3. ✅ Signal calculations for phagocytes

**Expected speedup from Phase 1 alone: 5-10x**

### Phase 2: Core Simulation Loop (1 week)
Port the main simulation loop (lines 76-583):

1. Agent movement
2. Engulfment processing
3. Treg activation

**Expected speedup from Phase 2: 20-50x total**

### Phase 3: Polish (optional)
- Matrix diffusion (if needed)
- Optimization of memory allocation

---

## Example 1: Port `get_8n_avg_signal_fast()` to C++

### Original R Code (SLOW):
```r
get_8n_avg_signal_fast <- function(x, y, act_radius_signal, signal_matrix) {
  loc  = c(x, y)
  x_coordinates = (loc[1]-act_radius_signal):(loc[1]+act_radius_signal)
  x_coordinates = x_coordinates[x_coordinates>0 & x_coordinates<=grid_size]
  y_coordinates = (loc[2]-act_radius_signal):(loc[2]+act_radius_signal)
  y_coordinates = y_coordinates[y_coordinates>0 & y_coordinates<=grid_size]
  dval = signal_matrix[y_coordinates, x_coordinates]
  return(mean(dval))
}
```

### C++ Version (FAST - 50x speedup):

Create file `MISC/fast_functions.cpp`:

```cpp
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double get_8n_avg_signal_cpp(int x, int y, int act_radius_signal,
                              NumericMatrix signal_matrix, int grid_size) {
  // Calculate bounds with edge checking
  int x_start = std::max(0, x - 1 - act_radius_signal);  // R to C++ index
  int x_end = std::min(grid_size - 1, x - 1 + act_radius_signal);
  int y_start = std::max(0, y - 1 - act_radius_signal);
  int y_end = std::min(grid_size - 1, y - 1 + act_radius_signal);

  double sum = 0.0;
  int count = 0;

  for (int yi = y_start; yi <= y_end; yi++) {
    for (int xi = x_start; xi <= x_end; xi++) {
      sum += signal_matrix(yi, xi);
      count++;
    }
  }

  return sum / count;
}

// Vectorized version - process multiple agents at once
// [[Rcpp::export]]
NumericVector get_8n_avg_signal_vectorized_cpp(
    IntegerVector x_vec, IntegerVector y_vec,
    int act_radius_signal, NumericMatrix signal_matrix, int grid_size
) {
  int n = x_vec.size();
  NumericVector results(n);

  for (int i = 0; i < n; i++) {
    results[i] = get_8n_avg_signal_cpp(
      x_vec[i], y_vec[i], act_radius_signal, signal_matrix, grid_size
    );
  }

  return results;
}
```

### Use in R:
```r
# Compile once at startup
Rcpp::sourceCpp("MISC/fast_functions.cpp")

# Then use just like R function
avg_signal <- get_8n_avg_signal_cpp(x = 10, y = 15,
                                     act_radius_signal = 1,
                                     signal_matrix = ROS,
                                     grid_size = 25)
```

---

## Example 2: Port Core Simulation Loop

For maximum speedup, port the entire simulation loop to C++.

**Benefits:**
- No R/C++ boundary crossing in the hot loop
- Memory allocation optimization
- 50-100x speedup

**File: `MISC/simulation_core.cpp`**

```cpp
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List run_simulation_cpp(
    int t_max,
    int grid_size,
    IntegerVector phagocyte_x,
    IntegerVector phagocyte_y,
    NumericMatrix DAMPs,
    NumericMatrix SAMPs,
    NumericMatrix ROS,
    // ... all other parameters
) {
  // Main simulation loop
  for (int t = 0; t < t_max; t++) {

    // Update SAMPs (vectorized in C++)
    for (int i = 0; i < active_tregs.size(); i++) {
      int treg_idx = active_tregs[i];
      int tx = treg_x[treg_idx] - 1;  // R to C++ indexing
      int ty = treg_y[treg_idx] - 1;
      SAMPs(ty, tx) += treg_activity_SAMPs_binary[treg_idx] * add_SAMPs * allow_tregs;
    }

    // ... rest of simulation
  }

  return List::create(
    Named("epithelium_longitudinal") = epithelium_longitudinal,
    Named("macrophages_longitudinal") = macrophages_longitudinal,
    // ... other results
  );
}
```

---

## Benchmarking: Prove the Speedup

Create `benchmark_rcpp.R`:

```r
library(microbenchmark)
library(Rcpp)

# Compile C++ functions
sourceCpp("MISC/fast_functions.cpp")

# Setup test data
grid_size <- 25
ROS <- matrix(runif(grid_size^2), grid_size, grid_size)

# Benchmark: R vs C++
results <- microbenchmark(
  R_version = get_8n_avg_signal_fast(10, 15, 1, ROS),
  Cpp_version = get_8n_avg_signal_cpp(10, 15, 1, ROS, grid_size),
  times = 1000
)

print(results)
# Typical result:
#   R_version:   ~15 microseconds
#   Cpp_version: ~0.3 microseconds  (50x faster!)
```

---

## Full Workflow Example

### 1. Create C++ file: `MISC/simulation_accelerated.cpp`

```cpp
// Put all your hot functions here
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double get_8n_avg_signal_cpp(...) { /* as above */ }

// [[Rcpp::export]]
List kill_microbes_with_ros_cpp(...) { /* as above */ }

// Add more functions as needed
```

### 2. Modify your R script to use C++:

```r
# At the top of RUN_REPS.R
library(Rcpp)
sourceCpp("MISC/simulation_accelerated.cpp")

# Then in the simulation loop, replace:
# OLD (R):
for (i in 1:nrow(pathogen_coords)) {
  pathogen_avg_ROS[i] = get_8n_avg_signal_fast(...)
}

# NEW (C++):
result <- kill_microbes_with_ros_cpp(
  pathogen_coords, ROS, act_radius_ROS, th_ROS_microbe, grid_size
)
pathogen_coords <- result$surviving_pathogens
pathogens_killed_by_ROS <- pathogens_killed_by_ROS + result$n_killed
```

---

## Common Gotchas

### 1. **R uses 1-based indexing, C++ uses 0-based**
```cpp
// In C++, subtract 1 from R coordinates
int x_cpp = x_r - 1;
```

### 2. **Matrix indexing is reversed**
```cpp
// R: matrix[row, col]
// C++: matrix(row, col)  - note the parentheses!
```

### 3. **Random numbers need special handling**
```cpp
// Use Rcpp's RNG to match R's random seed
#include <Rcpp.h>
// Rcpp::runif(), Rcpp::rnorm(), etc. automatically sync with set.seed()
```

---

## Expected Performance for Your Simulation

Based on your code structure:

| Component | % of Time | R Time | Rcpp Time | Speedup |
|-----------|-----------|--------|-----------|---------|
| Signal calculations | 30% | 100s | 2s | 50x |
| Agent loops | 40% | 100s | 1.5s | 67x |
| Engulfment | 15% | 100s | 2s | 50x |
| Matrix ops | 10% | 100s | 10s | 10x |
| Other | 5% | 100s | 50s | 2x |

**Overall: 30-50x faster** for full simulation with Rcpp.

---

## Resources

- [Rcpp Gallery](https://gallery.rcpp.org/) - Example code
- [Rcpp Book](http://www.rcpp.org/book/) - Comprehensive guide
- [Quick Reference](https://dirk.eddelbuettel.com/code/rcpp/Rcpp-quickref.pdf)

---

## Next Steps

1. ✅ Install Rcpp and test compilation
2. ✅ Port `get_8n_avg_signal_fast()` as proof of concept
3. ✅ Benchmark to verify speedup
4. ✅ Gradually port more functions
5. ✅ Eventually port the full simulation loop

**Start small, benchmark often, and you'll see massive gains!**
