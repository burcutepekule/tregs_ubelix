// Example: How to port the ROS microbe killing loop to C++
// This single function could give you 20-50x speedup on this operation

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List kill_microbes_with_ros_cpp(
    NumericMatrix pathogen_coords,
    NumericMatrix ROS,
    int act_radius_ROS,
    double th_ROS_microbe,
    int grid_size
) {
  int n_pathogens = pathogen_coords.nrow();
  std::vector<int> pathogens_to_keep;

  for (int i = 0; i < n_pathogens; i++) {
    int x = pathogen_coords(i, 0) - 1;  // R to C++ indexing
    int y = pathogen_coords(i, 1) - 1;

    // Calculate average ROS in neighborhood
    double sum_ros = 0.0;
    int count = 0;

    int x_start = std::max(0, x - act_radius_ROS);
    int x_end = std::min(grid_size - 1, x + act_radius_ROS);
    int y_start = std::max(0, y - act_radius_ROS);
    int y_end = std::min(grid_size - 1, y + act_radius_ROS);

    for (int yi = y_start; yi <= y_end; yi++) {
      for (int xi = x_start; xi <= x_end; xi++) {
        sum_ros += ROS(yi, xi);
        count++;
      }
    }

    double avg_ros = sum_ros / count;

    // Keep pathogen if it survives
    if (avg_ros <= th_ROS_microbe) {
      pathogens_to_keep.push_back(i);
    }
  }

  // Build result matrix
  int n_survivors = pathogens_to_keep.size();
  NumericMatrix survivors(n_survivors, 3);

  for (int i = 0; i < n_survivors; i++) {
    int idx = pathogens_to_keep[i];
    survivors(i, 0) = pathogen_coords(idx, 0);
    survivors(i, 1) = pathogen_coords(idx, 1);
    survivors(i, 2) = pathogen_coords(idx, 2);
  }

  colnames(survivors) = CharacterVector::create("x", "y", "id");

  return List::create(
    Named("surviving_pathogens") = survivors,
    Named("n_killed") = n_pathogens - n_survivors
  );
}

// Compile this with: Rcpp::sourceCpp("example_rcpp_speedup.cpp")
// Then call from R: result <- kill_microbes_with_ros_cpp(pathogen_coords, ROS, 1, 0.5, 25)
