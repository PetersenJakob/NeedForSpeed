#pragma once

#include <Eigen/Dense>
#include <vector>


namespace grid {

	std::vector<double> uniform(
		const double x_min, 
		const double x_max, 
		const int n_points);

	std::vector<double> exponential(
		const double x_min,
		const double x_max,
		const int n_points,
		const double scaling = 2.0);

	std::vector<double> hyperbolic(
		const double x_min,
		const double x_max,
		const int n_points,
		const double x_center = 0.0,
		const double scaling = 0.1);

	namespace eigen {

		Eigen::VectorXd uniform(
			const double x_min,
			const double x_max,
			const int n_points);

	}

}
