#include <functional>
#include <vector>

#include "convergence.h"


void grid_increment(
	const int n_increments,
	std::vector<double>& grid,
	std::function<std::vector<double>(double, double, int)> grid_generator) {

	grid = grid_generator(
		grid.front(), 
		grid.back(), 
		(int)grid.size() + n_increments);

}


std::vector<std::vector<double>> norm_vector(const int n_iterations) {

	std::vector<double> inner(n_iterations, 0.0);

	//
	// step_size_vec
	// max_norm
	// l1_vec_norm
	// l2_vec_norm
	// l1_func_norm
	// l2_func_norm
	//

	std::vector<std::vector<double>> outer(6, inner);

	return outer;

}


double average_grid_spacing(const std::vector<double>& grid) {

	return (grid.back() - grid.front()) / (grid.size() - 1);

}
