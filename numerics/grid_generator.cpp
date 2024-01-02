#include <vector>


// TODO: "Doc string" for function
// TODO: Description of parameters and return value type
std::vector<double> grid_equidistant(const double x_min, const double x_max, const int n_points) {

	// Constant step size.
	const double dx = (x_max - x_min) / (n_points - 1);

	//
	std::vector<double> grid(n_points);

	for (int n = 0; n != n_points; ++n)
		grid.at(n) = dx * n + x_min;

	return grid;
}
