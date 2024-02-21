#pragma once

#include <vector>


namespace test_util {

	std::vector<double> test_function(
		const std::vector<double>& grid,
		const int function_type,
		const int derivative_order);

	void print_test_results(
		const std::vector<double>& grid,
		const std::vector<double>& func,
		const std::vector<double>& deriv,
		const std::vector<double>& fd_result);

}
