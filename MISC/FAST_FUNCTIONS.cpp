// ============================================================================
// HIGH-PERFORMANCE C++ IMPLEMENTATIONS FOR TREGS SIMULATION
// ============================================================================
// This file contains C++ implementations of computational bottlenecks
// Expected speedup: 20-100x faster than R for these operations
//
// Compile with: Rcpp::sourceCpp("MISC/FAST_FUNCTIONS.cpp")
// ============================================================================

#include <Rcpp.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

// ============================================================================
// SIGNAL CALCULATIONS (MAJOR BOTTLENECK)
// ============================================================================

// [[Rcpp::export]]
double get_8n_avg_signal_cpp(int x, int y, int act_radius_signal,
                             NumericMatrix signal_matrix, int grid_size)
{
  // Convert from R's 1-based to C++'s 0-based indexing
  int x0 = x - 1;
  int y0 = y - 1;

  // Calculate bounds with edge checking
  int x_start = std::max(0, x0 - act_radius_signal);
  int x_end = std::min(grid_size - 1, x0 + act_radius_signal);
  int y_start = std::max(0, y0 - act_radius_signal);
  int y_end = std::min(grid_size - 1, y0 + act_radius_signal);

  double sum = 0.0;
  int count = 0;

  for (int yi = y_start; yi <= y_end; yi++)
  {
    for (int xi = x_start; xi <= x_end; xi++)
    {
      sum += signal_matrix(yi, xi);
      count++;
    }
  }

  return (count > 0) ? (sum / count) : 0.0;
}

// Vectorized version - process multiple agents at once
// [[Rcpp::export]]
NumericVector get_8n_avg_signal_vectorized_cpp(
    IntegerVector x_vec, IntegerVector y_vec,
    int act_radius_signal, NumericMatrix signal_matrix, int grid_size)
{
  int n = x_vec.size();
  NumericVector results(n);

  for (int i = 0; i < n; i++)
  {
    results[i] = get_8n_avg_signal_cpp(
        x_vec[i], y_vec[i], act_radius_signal, signal_matrix, grid_size);
  }

  return results;
}

// ============================================================================
// MICROBE KILLING WITH ROS (VERY HOT LOOP)
// ============================================================================

// [[Rcpp::export]]
List kill_microbes_with_ros_cpp(
    NumericMatrix microbe_coords,
    NumericMatrix ROS,
    int act_radius_ROS,
    double th_ROS_microbe,
    int grid_size)
{
  int n_microbes = microbe_coords.nrow();

  if (n_microbes == 0)
  {
    NumericMatrix empty(0, 3);
    colnames(empty) = CharacterVector::create("x", "y", "id");
    return List::create(
        Named("surviving_microbes") = empty,
        Named("n_killed") = 0);
  }

  std::vector<int> survivors;
  survivors.reserve(n_microbes); // Pre-allocate for efficiency

  for (int i = 0; i < n_microbes; i++)
  {
    int x = static_cast<int>(microbe_coords(i, 0));
    int y = static_cast<int>(microbe_coords(i, 1));

    // Calculate average ROS using the fast function
    double avg_ros = get_8n_avg_signal_cpp(x, y, act_radius_ROS, ROS, grid_size);

    // Keep microbe if it survives
    if (avg_ros <= th_ROS_microbe)
    {
      survivors.push_back(i);
    }
  }

  // Build result matrix with survivors
  int n_survivors = survivors.size();
  NumericMatrix result(n_survivors, 3);

  for (int i = 0; i < n_survivors; i++)
  {
    int idx = survivors[i];
    result(i, 0) = microbe_coords(idx, 0);
    result(i, 1) = microbe_coords(idx, 1);
    result(i, 2) = microbe_coords(idx, 2);
  }

  colnames(result) = CharacterVector::create("x", "y", "id");

  return List::create(
      Named("surviving_microbes") = result,
      Named("n_killed") = n_microbes - n_survivors);
}

// ============================================================================
// PHAGOCYTE SIGNAL CALCULATIONS (HOT LOOP)
// ============================================================================

// Calculate DAMPs and SAMPs for all specified phagocytes
// [[Rcpp::export]]
List calculate_phagocyte_signals_cpp(
    IntegerVector phagocyte_indices,
    IntegerVector phagocyte_x,
    IntegerVector phagocyte_y,
    NumericMatrix phagocyte_bacteria_registry,
    int act_radius_DAMPs,
    int act_radius_SAMPs,
    NumericMatrix DAMPs,
    NumericMatrix SAMPs,
    int grid_size)
{
  int n = phagocyte_indices.size();
  NumericVector avg_DAMPs(n);
  NumericVector avg_SAMPs(n);
  IntegerVector bacteria_counts(n);

  for (int i = 0; i < n; i++)
  {
    int idx = phagocyte_indices[i] - 1; // R to C++ indexing

    // Get position
    int x = phagocyte_x[idx];
    int y = phagocyte_y[idx];

    // Calculate average signals
    avg_DAMPs[i] = get_8n_avg_signal_cpp(x, y, act_radius_DAMPs, DAMPs, grid_size);
    avg_SAMPs[i] = get_8n_avg_signal_cpp(x, y, act_radius_SAMPs, SAMPs, grid_size);

    // Count bacteria
    int bacteria_count = 0;
    for (int j = 0; j < phagocyte_bacteria_registry.ncol(); j++)
    {
      bacteria_count += phagocyte_bacteria_registry(idx, j);
    }
    bacteria_counts[i] = bacteria_count;
  }

  return List::create(
      Named("avg_DAMPs") = avg_DAMPs,
      Named("avg_SAMPs") = avg_SAMPs,
      Named("bacteria_counts") = bacteria_counts);
}

// ============================================================================
// ENGULFMENT PROCESSING (COMPLEX BOTTLENECK)
// ============================================================================

// Process engulfment for a single phagocyte
// Returns indices of microbes to remove
// [[Rcpp::export]]
IntegerVector find_overlapping_microbes_cpp(
    int px, int py,
    NumericMatrix microbe_coords)
{
  int n_microbes = microbe_coords.nrow();
  std::vector<int> overlapping;

  for (int i = 0; i < n_microbes; i++)
  {
    if (microbe_coords(i, 0) == px && microbe_coords(i, 1) == py)
    {
      overlapping.push_back(i + 1); // Convert to R indexing
    }
  }

  return wrap(overlapping);
}

// ============================================================================
// TREG VICINITY CALCULATIONS (HOT LOOP)
// ============================================================================

// Find tregs within vicinity of a phagocyte
// [[Rcpp::export]]
IntegerVector find_nearby_tregs_cpp(
    int px, int py,
    IntegerVector treg_x,
    IntegerVector treg_y,
    int treg_vicinity_effect)
{
  int n_tregs = treg_x.size();
  std::vector<int> nearby_indices;

  for (int i = 0; i < n_tregs; i++)
  {
    int dist_x = std::abs(treg_x[i] - px);
    int dist_y = std::abs(treg_y[i] - py);

    if (dist_x <= treg_vicinity_effect && dist_y <= treg_vicinity_effect)
    {
      nearby_indices.push_back(i + 1); // Convert to R indexing
    }
  }

  return wrap(nearby_indices);
}

// ============================================================================
// EPITHELIAL ROS CALCULATIONS (VECTORIZABLE)
// ============================================================================

// Calculate mean ROS for all epithelial cells at once
// [[Rcpp::export]]
NumericVector calculate_epithelial_ros_cpp(
    IntegerVector epithelium_x,
    int act_radius_ROS,
    NumericMatrix ROS,
    int grid_size)
{
  int n_cells = epithelium_x.size();
  NumericVector ros_means(n_cells);

  for (int i = 0; i < n_cells; i++)
  {
    int px = epithelium_x[i];
    int x0 = px - 1; // Convert to 0-based indexing

    // Calculate x coordinate range
    int x_start = std::max(0, x0 - act_radius_ROS);
    int x_end = std::min(grid_size - 1, x0 + act_radius_ROS);

    // ROS is at row 0 (epithelium is at y=1 in R, index 0 in C++)
    double sum = 0.0;
    int count = 0;

    for (int xi = x_start; xi <= x_end; xi++)
    {
      sum += ROS(0, xi); // Row 0 is the epithelium layer
      count++;
    }

    ros_means[i] = (count > 0) ? (sum / count) : 0.0;
  }

  return ros_means;
}

// ============================================================================
// MATRIX DIFFUSION (IF NEEDED)
// ============================================================================

// Optimized 8-neighbor diffusion
// [[Rcpp::export]]
NumericMatrix diffuse_matrix_cpp(
    NumericMatrix mat,
    double D,
    double max_cell_value)
{
  int nr = mat.nrow();
  int nc = mat.ncol();

  // Create padded matrix
  NumericMatrix padded(nr + 2, nc + 2);
  for (int i = 0; i < nr; i++)
  {
    for (int j = 0; j < nc; j++)
    {
      padded(i + 1, j + 1) = mat(i, j);
    }
  }

  // Calculate 8-neighbor Laplacian
  NumericMatrix result(nr, nc);

  for (int i = 0; i < nr; i++)
  {
    for (int j = 0; j < nc; j++)
    {
      double laplacian =
          padded(i, j) +         // top-left
          padded(i, j + 1) +     // top
          padded(i, j + 2) +     // top-right
          padded(i + 1, j) +     // left
          padded(i + 1, j + 2) + // right
          padded(i + 2, j) +     // bottom-left
          padded(i + 2, j + 1) + // bottom
          padded(i + 2, j + 2) - // bottom-right
          8.0 * mat(i, j);       // center

      double new_val = mat(i, j) + D * laplacian;
      result(i, j) = std::min(max_cell_value, new_val);
    }
  }

  return result;
}

// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================

// Fast shift and insert for bacteria registry
// [[Rcpp::export]]
NumericVector shift_insert_fast_cpp(
    NumericVector vec,
    NumericVector insert_vals)
{
  int n_insert = insert_vals.size();
  int n_vec = vec.size();

  NumericVector result(n_vec);

  if (n_insert >= n_vec)
  {
    // If inserting more than vec size, just take first n_vec elements
    for (int i = 0; i < n_vec; i++)
    {
      result[i] = insert_vals[i];
    }
  }
  else
  {
    // Insert at beginning, shift rest
    for (int i = 0; i < n_insert; i++)
    {
      result[i] = insert_vals[i];
    }
    for (int i = n_insert; i < n_vec; i++)
    {
      result[i] = vec[i - n_insert];
    }
  }

  return result;
}

// Optimized coordinate movement helper
// [[Rcpp::export]]
IntegerVector iszero_coordinates_cpp(IntegerVector x)
{
  int n = x.size();
  IntegerVector y(n);

  // First pass: sample for all positions
  NumericVector rand_vals = runif(n, -1, 1);

  for (int i = 0; i < n; i++)
  {
    if (x[i] == 0)
    {
      // For zero values, choose -1 or 1
      y[i] = (rand_vals[i] < 0) ? -1 : 1;
    }
    else
    {
      // For non-zero, sample from {-1, 0, 1}
      double r = runif(1)[0];
      if (r < 0.333)
      {
        y[i] = -1;
      }
      else if (r < 0.667)
      {
        y[i] = 0;
      }
      else
      {
        y[i] = 1;
      }
    }
  }

  return y;
}

// ============================================================================
// BATCH OPERATIONS FOR MAXIMUM EFFICIENCY
// ============================================================================

// Update SAMPs matrix with multiple active tregs at once
// [[Rcpp::export]]
NumericMatrix update_SAMPs_batch_cpp(
    NumericMatrix SAMPs,
    IntegerVector active_tregs,
    IntegerVector treg_x,
    IntegerVector treg_y,
    NumericVector treg_activity_SAMPs_binary,
    double add_SAMPs,
    int allow_tregs)
{
  NumericMatrix result = clone(SAMPs);
  int n_active = active_tregs.size();

  for (int i = 0; i < n_active; i++)
  {
    int idx = active_tregs[i] - 1; // R to C++ indexing
    int x = treg_x[idx] - 1;
    int y = treg_y[idx] - 1;

    result(y, x) += treg_activity_SAMPs_binary[idx] * add_SAMPs * allow_tregs;
  }

  return result;
}
