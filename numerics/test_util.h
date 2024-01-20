#pragma once

#include <vector>

namespace test_util {

	std::vector<double> test_function(
		std::vector<double> grid,
		const int function_type,
		const int derivative_order);


	void print_test_results(
		std::vector<double> grid,
		std::vector<double> func,
		std::vector<double> deriv,
		std::vector<double> fd_result);


	// Element-wise subtraction of vectors.
	std::vector<double> vector_diff(
		const std::vector<double>& vec1,
		const std::vector<double>& vec2);


	// Max norm.
	double max_norm(std::vector<double> vec);


	// l2-norm (vector norm).
	double l2_norm(const double dx, std::vector<double> vec);

}