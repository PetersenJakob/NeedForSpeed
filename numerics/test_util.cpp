#define _USE_MATH_DEFINES
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>

#include "test_util.h"


std::vector<double> test_util::test_function(
	std::vector<double> grid,
	const int function_type,
	const int derivative_order) {

	std::vector<double> inner(grid.size(), 0.0);

	std::vector<std::vector<double>> f_exp(5, inner);
	std::vector<std::vector<double>> f_cos(5, inner);
	std::vector<std::vector<double>> f_sum(5, inner);

	for (int i = 0; i != grid.size(); ++i) {

		f_exp[4][i] = 4.0 * exp(2.0 * grid[i]);

		f_exp[3][i] = 2.0 * exp(2.0 * grid[i]);

		f_exp[2][i] = exp(2.0 * grid[i]);

		f_exp[1][i] = exp(2.0 * grid[i]) / 2.0;

		f_exp[0][i] = exp(2.0 * grid[i]) / 4.0;


		f_cos[4][i] = -pow(M_PI, 2.0) * cos(M_PI * grid[i]);

		f_cos[3][i] = -M_PI * sin(M_PI * grid[i]);

		f_cos[2][i] = cos(M_PI * grid[i]);

		f_cos[1][i] = sin(M_PI * grid[i]) / M_PI;

		f_cos[0][i] = -cos(M_PI * grid[i]) / pow(M_PI, 2.0);


		f_sum[4][i] = f_cos[4][i] + f_exp[4][i];

		f_sum[3][i] = f_cos[3][i] + f_exp[3][i];

		f_sum[2][i] = f_cos[2][i] + f_exp[2][i];

		f_sum[1][i] = f_cos[1][i] + f_exp[1][i];

		f_sum[0][i] = f_cos[0][i] + f_exp[0][i];

	}

	const int function_idx = derivative_order + 2;

	if (function_type == 0) {
		return f_exp[function_idx];
	}
	else if (function_type == 1) {
		return f_cos[function_idx];
	}
	else {
		return f_sum[function_idx];
	}

}


void test_util::print_test_results(
	std::vector<double> grid,
	std::vector<double> func,
	std::vector<double> deriv,
	std::vector<double> fd_result) {

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


// l2-norm (vector norm). TODO: Should this be dependent on dx? Like a "Riemann sum" expression?
double test_util::l2_norm(const double dx, std::vector<double> vec) {

	// Square of each element.
	std::transform(vec.begin(), vec.end(), vec.begin(), [dx](double x) { return dx * pow(x, 2); });

	return sqrt(std::accumulate(vec.begin(), vec.end(), (double)0));

}
