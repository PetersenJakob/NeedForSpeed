#pragma once

#include <Eigen/Dense>
#include <vector>


// References:
// - Veldman and Rinzeman (1992)
// - Bodeau, Riboulet and Roncalli (2000)
// - White (2013)

// Standard form of 1-dimensional grid function:
// std::vector<double> grid::func(
//     const double x_min, 
//     const double x_max, 
//     const int n_points)
namespace grid {

	std::vector<double> uniform(
		const double x_min, 
		const double x_max, 
		const int n_points);

	std::vector<double> exponential_full(
		const double x_min,
		const double x_max,
		const int n_points,
		const double scaling = 2.0);

	std::vector<double> exponential(
		const double x_min,
		const double x_max,
		const int n_points);

	std::vector<double> hyperbolic_full(
		const double x_min,
		const double x_max,
		const int n_points,
		const double x_center = 0.0,
		const double scaling = 0.1);

	std::vector<double> hyperbolic(
		const double x_min,
		const double x_max,
		const int n_points);

	namespace eigen {

		Eigen::VectorXd uniform(
			const double x_min,
			const double x_max,
			const int n_points);

	}

}
