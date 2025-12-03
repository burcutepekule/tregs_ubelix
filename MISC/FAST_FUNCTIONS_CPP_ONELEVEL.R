# ============================================================================
# C++ ACCELERATED FUNCTIONS FOR TREGS SIMULATION
# ============================================================================
# This file loads the C++ implementations and provides R wrappers
# Expected speedup: 20-100x faster than pure R
# ============================================================================

library(Rcpp)

# ============================================================================
# COMPILE AND LOAD C++ FUNCTIONS
# ============================================================================

# Get the directory where this script is located
script_dir <- dirname(sys.frame(1)$ofile)
if (length(script_dir) == 0 || script_dir == "") {
  script_dir <- getwd()
}

cpp_file <- file.path(script_dir, "FAST_FUNCTIONS_ONELEVEL.cpp")

# Check if C++ file exists
if (!file.exists(cpp_file)) {
  stop(paste0("C++ file not found: ", cpp_file,
              "\nPlease ensure FAST_FUNCTIONS.cpp is in the MISC directory"))
}

# Compile and load C++ functions
cat("Compiling C++ functions... ")
tryCatch({
  sourceCpp(cpp_file)
  cat("✓ SUCCESS\n")
  cat("C++ acceleration enabled (20-100x speedup expected)\n")
  USE_CPP <- TRUE
}, error = function(e) {
  cat("✗ FAILED\n")
  cat("Error:", conditionMessage(e), "\n")
  cat("Falling back to R implementations (slower)\n")
  USE_CPP <- FALSE
})

# ============================================================================
# ORIGINAL R FUNCTIONS (FALLBACK IF C++ FAILS)
# ============================================================================

# Logistic function
logistic_function <- function(x, k = k_in, x0 = x0_in) {
  return(1 / (1 + exp(-k * (x - x0))))
}

logistic_scaled_0_to_5_quantized <- function(x,  k = k_in, x0 = x0_in) {
  return(round(5*plogis(x, location = x0, scale = 1 / k)))
}

# Original R implementation (kept as fallback)
get_8n_avg_signal_fast <- function(x, y, act_radius_signal, signal_matrix) {
  loc  = c(x, y)
  x_coordinates = (loc[1]-act_radius_signal):(loc[1]+act_radius_signal)
  x_coordinates = x_coordinates[x_coordinates>0 & x_coordinates<=grid_size]
  y_coordinates = (loc[2]-act_radius_signal):(loc[2]+act_radius_signal)
  y_coordinates = y_coordinates[y_coordinates>0 & y_coordinates<=grid_size]
  dval = signal_matrix[y_coordinates, x_coordinates]
  return(mean(dval))
}

# Fast shift_insert for matrix operations
shift_insert_fast <- function(vec, insert_vals) {
  n_insert <- length(insert_vals)
  n_vec <- length(vec)
  
  if(n_insert >= n_vec) {
    return(insert_vals[1:n_vec])
  } else {
    return(c(insert_vals, vec[1:(n_vec - n_insert)]))
  }
}

# Optimized version of iszero_coordinates
iszero_coordinates <- function(x) {
  y <- sample(c(-1, 0, 1), length(x), replace = TRUE)
  zero_idx <- which(x == 0)
  y[zero_idx] <- sample(c(-1, 1), length(zero_idx), replace = TRUE)
  return(y)
}

# Diffusion matrix (R version)
diffuse_matrix <- function(mat, D, max_cell_value) {
  nr <- nrow(mat)
  nc <- ncol(mat)
  
  padded <- matrix(0, nrow = nr + 2, ncol = nc + 2)
  padded[2:(nr + 1), 2:(nc + 1)] <- mat
  
  laplacian <- (
    padded[1:nr,     1:nc    ] +
      padded[1:nr,     2:(nc+1)] +
      padded[1:nr,     3:(nc+2)] +
      padded[2:(nr+1), 1:nc    ] +
      padded[2:(nr+1), 3:(nc+2)] +
      padded[3:(nr+2), 1:nc    ] +
      padded[3:(nr+2), 2:(nc+1)] +
      padded[3:(nr+2), 3:(nc+2)]
    - 8 * mat
  )
  
  mat_new <- mat + D * laplacian
  mat_new <- matrix(pmin(max_cell_value, mat_new), nrow = nrow(mat), ncol = ncol(mat))
  
  return(mat_new)
}

shift_insert <- function(current_registry, new_elements_vector) {
  combined_registry <- c(new_elements_vector, current_registry)
  target_length <- length(current_registry)
  result_registry <- combined_registry[1:target_length]
  
  return(result_registry)
}

get_middle_percent <- function(seq_vector, percent) {
  n_total <- length(seq_vector)
  n_select <- ceiling(n_total * percent / 100)
  
  mid <- floor(n_total / 2)
  half_window <- floor(n_select / 2)
  
  start_idx <- max(1, mid - half_window + 1)
  end_idx <- min(n_total, start_idx + n_select - 1)
  
  return(seq_vector[start_idx:end_idx])
}

sample_rbeta <- function(alpha, beta) {
  x <- rgamma(1, shape = alpha, rate = 1.0)
  y <- rgamma(1, shape = beta, rate = 1.0)
  return(x / (x + y))
}

# ============================================================================
# WRAPPER FUNCTIONS (USE C++ IF AVAILABLE, OTHERWISE R)
# ============================================================================

# Wrapper for get_8n_avg_signal - uses C++ if available
get_8n_avg_signal_fast_wrapper <- function(x, y, act_radius_signal, signal_matrix) {
  if (exists("get_8n_avg_signal_cpp", mode = "function")) {
    return(get_8n_avg_signal_cpp(x, y, act_radius_signal, signal_matrix, grid_size))
  } else {
    return(get_8n_avg_signal_fast(x, y, act_radius_signal, signal_matrix))
  }
}

# Vectorized version
get_8n_avg_signal_vectorized <- function(x_vec, y_vec, act_radius_signal, signal_matrix, grid_size) {
  if (exists("get_8n_avg_signal_vectorized_cpp", mode = "function")) {
    return(get_8n_avg_signal_vectorized_cpp(x_vec, y_vec, act_radius_signal, signal_matrix, grid_size))
  } else {
    # Fallback to R implementation
    n_agents = length(x_vec)
    results = numeric(n_agents)
    
    for (i in 1:n_agents) {
      x = x_vec[i]
      y = y_vec[i]
      
      x_coordinates = (x - act_radius_signal):(x + act_radius_signal)
      x_coordinates = x_coordinates[x_coordinates > 0 & x_coordinates <= grid_size]
      
      y_coordinates = (y - act_radius_signal):(y + act_radius_signal)
      y_coordinates = y_coordinates[y_coordinates > 0 & y_coordinates <= grid_size]
      
      dval = signal_matrix[y_coordinates, x_coordinates]
      results[i] = mean(dval)
    }
    
    return(results)
  }
}

# ============================================================================
# INFORMATION FUNCTIONS
# ============================================================================

# Check which functions are using C++
check_cpp_status <- function() {
  cpp_functions <- c(
    "get_8n_avg_signal_cpp",
    "get_8n_avg_signal_vectorized_cpp",
    "kill_microbes_with_ros_cpp",
    "calculate_phagocyte_signals_cpp",
    "find_overlapping_microbes_cpp",
    "find_nearby_tregs_cpp",
    "calculate_epithelial_ros_cpp",
    "diffuse_matrix_cpp",
    "shift_insert_fast_cpp",
    "iszero_coordinates_cpp",
    "update_SAMPs_batch_cpp",
    "update_ROS_batch_cpp"
  )
  
  cat("\n")
  cat("=" = rep("=", 70), sep = "")
  cat("\nC++ ACCELERATION STATUS\n")
  cat("=" = rep("=", 70), sep = "")
  cat("\n\n")
  
  n_available <- 0
  for (func_name in cpp_functions) {
    if (exists(func_name, mode = "function")) {
      cat("✓", func_name, "\n")
      n_available <- n_available + 1
    } else {
      cat("✗", func_name, "(using R fallback)\n")
    }
  }
  
  cat("\n")
  cat("Total C++ functions available:", n_available, "/", length(cpp_functions), "\n")
  
  if (n_available == length(cpp_functions)) {
    cat("\n🚀 Full C++ acceleration active! (Expected 20-100x speedup)\n")
  } else if (n_available > 0) {
    cat("\n⚠ Partial C++ acceleration (some functions using R fallback)\n")
  } else {
    cat("\n⚠ No C++ acceleration (all functions using R - will be slower)\n")
    cat("   Check compilation errors above\n")
  }
  
  cat("=" = rep("=", 70), sep = "")
  cat("\n\n")
}

# Print status when file is loaded
cat("\n")
cat("FAST_FUNCTIONS_CPP.R loaded successfully\n")
cat("Use check_cpp_status() to see which functions are accelerated\n")
cat("\n")
