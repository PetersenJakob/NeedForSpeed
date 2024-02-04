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


// Element-wise subtraction of vectors.
std::vector<double> test_util::vector_diff(
	const std::vector<double>& vec1,
	const std::vector<double>& vec2) {

	std::vector<double> diff(vec1.size(), 0.0);

	std::transform(vec1.begin(), vec1.end(), vec2.begin(), diff.begin(), std::minus<double>());

	return diff;

}


// Max norm.
double test_util::max_norm(std::vector<double> vec) {

	// Absolue value of each element.
	std::transform(vec.begin(), vec.end(), vec.begin(), [](double x) { return abs(x); });

	return *std::max_element(vec.begin(), vec.end());

}


// l2-norm (vector norm). 
// TODO: Should this be dependent on dx? Like a "Riemann sum" expression?
double test_util::l2_norm(const double dx, std::vector<double> vec) {

	// Square of each element.
//	std::transform(vec.begin(), vec.end(), vec.begin(), [dx](double x) { return dx * x * x; });

//	return sqrt(std::accumulate(vec.begin(), vec.end(), (double)0));

	double result = 0.0;

	for (int i = 0; i != vec.size() - 1; ++i) {

		result += dx * abs(vec[i + 1] + vec[i]) / 2.0;
//		result += abs(vec[i + 1] + vec[i]) / 2.0;
//		result += vec[i];

	}

	return result;

}


// l2-norm (vector norm). 
// TODO: Should this be dependent on dx? Like a "Riemann sum" expression?
double test_util::l2_norm(const std::vector<double>& grid, std::vector<double> vec) {

	// Square of each element.
	std::transform(vec.begin(), vec.end(), grid.begin(), vec.begin(), [](double x, double dx) { return dx * x * x; });

	return sqrt(std::accumulate(vec.begin(), vec.end(), (double)0));

}


// Simple linear regression.
std::vector<double> test_util::slr(const std::vector<double>& x, const std::vector<double>& y) {

	double x_mean = 0.0;
	double y_mean = 0.0;

	double xx = 0.0;
	double xy = 0.0;

	for (int i = 0; i != x.size(); ++i) {

		x_mean += x[i];
		y_mean += y[i];

	}

	x_mean /= x.size();
	y_mean /= y.size();

	for (int i = 0; i != x.size(); ++i) {

		xx += pow(x[i] - x_mean, 2);
		xy += (x[i] - x_mean) * (y[i] - y_mean);

	}

	double beta = xy / xx;
	double alpha = y_mean - beta * x_mean;

	std::vector<double> result{ beta, alpha };

	return result;

}
