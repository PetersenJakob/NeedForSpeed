#define _USE_MATH_DEFINES
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>

#include "test_util.h"


std::vector<double> test_util::test_function(
	const std::vector<double>& grid,
	const int function_type,
	const int derivative_order) {

	// Function type:
	//	 0: Exponential function.
	//	 1: Cosine function.
	//	 2: Sum of exponential and cosine.

	// Derivative order:
	//	 2: 2nd order derivative.
	//	 1: 1st order derivative.
	//	 0: Function.
	//	-1: 1st order antiderivative.
	//	-2: 2nd order antiderivative.

	const int function_idx = derivative_order + 2;

	std::vector<double> inner(grid.size(), 0.0);

	std::vector<std::vector<double>> f_exp(5, inner);
	std::vector<std::vector<double>> f_cos(5, inner);
	std::vector<std::vector<double>> f_sum(5, inner);

	const double pi = M_PI;
	const double pi_sq = pi * pi;

	for (int i = 0; i != grid.size(); ++i) {

		// #####################
		// Exponential function.
		// #####################
		
		// 2nd order derivative.
		f_exp[4][i] = 4.0 * exp(2.0 * grid[i]);
		// 1st order derivative.
		f_exp[3][i] = 2.0 * exp(2.0 * grid[i]);
		// Function.
		f_exp[2][i] = exp(2.0 * grid[i]);
		// 1st order antiderivative.
		f_exp[1][i] = exp(2.0 * grid[i]) / 2.0;
		// 2nd order antiderivative.
		f_exp[0][i] = exp(2.0 * grid[i]) / 4.0;

		// ################
		// Cosine function.
		// ################

		// 2nd order derivative.
		f_cos[4][i] = -pi_sq * cos(pi * grid[i]);
		// 1st order derivative.
		f_cos[3][i] = -pi * sin(pi * grid[i]);
		// Function.
		f_cos[2][i] = cos(pi * grid[i]);
		// 1st order antiderivative.
		f_cos[1][i] = sin(pi * grid[i]) / pi;
		// 2nd order antiderivative.
		f_cos[0][i] = -cos(pi * grid[i]) / pi_sq;

		// ##############################
		// Sum of exponential and cosine.
		// ##############################

		// 2nd order derivative.
		f_sum[4][i] = f_exp[4][i] + f_cos[4][i];
		// 1st order derivative.
		f_sum[3][i] = f_exp[3][i] + f_cos[3][i];
		// Function.
		f_sum[2][i] = f_exp[2][i] + f_cos[2][i];
		// 1st order antiderivative.
		f_sum[1][i] = f_exp[1][i] + f_cos[1][i];
		// 2nd order antiderivative.
		f_sum[0][i] = f_exp[0][i] + f_cos[0][i];

	}

	if (function_type == 0) {
		return f_exp[function_idx];
	}
	else if (function_type == 1) {
		return f_cos[function_idx];
	}
	else if (function_type == 2) {
		return f_sum[function_idx];
	}
	else {
		throw std::invalid_argument("Function type unknown.");
	}

}


void test_util::print_test_results(
	const std::vector<double>& grid,
	const std::vector<double>& func,
	const std::vector<double>& deriv,
	const std::vector<double>& fd_result) {

	std::cout << std::scientific << std::setprecision(5);

	std::cout << std::endl 
		<< "             Grid          Func         Deriv      FD deriv      Abs diff" << std::endl;

	for (int i = 0; i != grid.size(); ++i) {

		std::cout
			<< std::setw(3) << i
			<< std::setw(14) << grid[i]
			<< std::setw(14) << func[i]
			<< std::setw(14) << deriv[i]
			<< std::setw(14) << fd_result[i]
			<< std::setw(14) << abs(deriv[i] - fd_result[i])
			<< std::endl;

	}

	std::cout << std::endl;

}
