#define _USE_MATH_DEFINES
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

#include "band_diagonal_matrix.h"
#include "finite_difference.h"
#include "grid_generator.h"
#include "tridiagonal_matrix_solver.h"


std::vector<double> vector_diff(
	const std::vector<double>& vec1,
	const std::vector<double>& vec2);

double max_norm(std::vector<double> vec);

double l2_norm(const double dx, std::vector<double> vec);


std::vector<double> test_exp(
	const std::vector<double>& grid, 
	const int deriv_order, 
	BandDiagonal& deriv_operator, 
	bool screen = false);

void test_cos(
	const std::vector<double>& grid, 
	const int deriv_order, 
	BandDiagonal& deriv_operator);

void print_test(
	const std::vector<double>& grid, 
	const std::vector<double>& func,
	const std::vector<double>& deriv,
	const std::vector<double>& deriv_fd);


int main() {

	const int order = 21;

	std::vector<double> grid_exp = grid_equidistant(-0.2, 0.2, order);
	std::vector<double> grid_cos = grid_equidistant(-M_PI / 2.0 , M_PI / 2.0, order);

	const double dx_exp = grid_exp[1] - grid_exp[0];
	const double dx_cos = grid_cos[1] - grid_cos[0];


	TriDiagonal d1_c2b1_exp = d1dx1::c2b1(order, dx_exp);
	TriDiagonal d1_c2b1_cos = d1dx1::c2b1(order, dx_cos);

	TriDiagonal d1_c2b2_exp = d1dx1::c2b2(order, dx_exp);
	TriDiagonal d1_c2b2_cos = d1dx1::c2b2(order, dx_cos);

	PentaDiagonal d1_c4b4_exp = d1dx1::c4b4(order, dx_exp);
	PentaDiagonal d1_c4b4_cos = d1dx1::c4b4(order, dx_cos);

	TriDiagonal d2_c2b1_exp = d2dx2::c2b1(order, dx_exp);
	TriDiagonal d2_c2b1_cos = d2dx2::c2b1(order, dx_cos);

	TriDiagonal d2_c2b2_exp = d2dx2::c2b2(order, dx_exp);
	TriDiagonal d2_c2b2_cos = d2dx2::c2b2(order, dx_cos);

	PentaDiagonal d2_c4b4_exp = d2dx2::c4b4(order, dx_exp);
	PentaDiagonal d2_c4b4_cos = d2dx2::c4b4(order, dx_cos);


	TriDiagonal d1_c2b1_exp_non = d1dx1::nonequidistant::c2b1(order, grid_exp);
	TriDiagonal d1_c2b1_cos_non = d1dx1::nonequidistant::c2b1(order, grid_cos);

	print_matrix(d1_c2b1_cos_non);

	print_matrix(d2_c2b1_cos);


//	test_exp(grid_exp, 1, d1_c4b4_exp);

//	test_cos(grid_cos, 1, d1_c2b2_cos);

//	test_cos(grid_cos, 1, d1_c4b4_cos);

#if false
	for (int i = 0; i != 19; ++i) {
	
		const int n_points = 20 + 10 * i;
		const double x_min = -0.2;
		const double x_max = 0.2;

		std::vector<double> grid = grid_equidistant(x_min, x_max, n_points);

		const double dx = grid[1] - grid[0];

//		TriDiagonal deriv_operator = d1dx1::c2b2(n_points, dx);
		PentaDiagonal deriv_operator = d1dx1::c4b4(n_points, dx);

		test_exp(grid, 1, deriv_operator, true);
	
	}
#endif

	std::vector<double> column_tmp(order, 0.0);

	std::vector<double> func(order, 0.0);

	for (int i = 0; i != order; ++i) {
		func[i] = cos(grid_cos[i]);
	}

	std::vector<double> solution_fd = func;

	const double dt = 0.01;
	const int n_steps = 3;

	print_matrix(d2_c2b1_cos);

	for (int i = 0; i != d2_c2b1_cos.order(); ++i) {

		for (int j = 0; j != d2_c2b1_cos.n_diagonals(); ++j) {
			d2_c2b1_cos.matrix[j][i] *= -dt;
		}

		d2_c2b1_cos.matrix[d2_c2b1_cos.bandwidth()][i] += 1.0;

	}

	for (int i = 0; i != 2 * d2_c2b1_cos.n_boundary_rows(); ++i) {
		for (int j = 0; j != d2_c2b1_cos.n_boundary_elements(); ++j) {
			d2_c2b1_cos.boundary_rows[i][j] *= -dt;
		}
	}

	d2_c2b1_cos.boundary_rows[0][0] += 1.0;
	d2_c2b1_cos.boundary_rows[1][d2_c2b1_cos.n_boundary_elements() - 1] += 1.0;


	print_matrix(d2_c2b1_cos);
	d2_c2b1_cos.adjust_boundary(column_tmp);
	print_matrix(d2_c2b1_cos);

	std::cout << "Column vector: " << std::endl;
	for (int i = 0; i != column_tmp.size(); ++i) {
		std::cout << column_tmp[i] << std::endl;
	}
	std::cout << std::endl;


	for (int i = 0; i != n_steps; ++i) {
		tri_solver(d2_c2b1_cos, solution_fd);
	}


	std::cout << std::scientific << std::setprecision(5);

	
	for (int i = 0; i != order; ++i) {

		std::cout
			<< std::setw(3) << i
			<< std::setw(14) << grid_cos[i]
			<< std::setw(14) << func[i]
			<< std::setw(14) << exp(-n_steps * dt) * func[i]
			<< std::setw(14) << solution_fd[i]
			<< std::setw(14) << abs(exp(-n_steps * dt) * func[i] - solution_fd[i])
			<< std::endl;

	}

	return 0;

}


std::vector<double> test_exp(const std::vector<double>& grid, const int deriv_order, BandDiagonal& deriv_operator, bool screen) {

	const int order = grid.size();

	std::vector<double> func(order, 0.0);
	std::vector<double> deriv1(order, 0.0);
	std::vector<double> deriv2(order, 0.0);

	for (int i = 0; i != order; ++i) {
		func[i] = exp(2 * grid[i]);
		deriv1[i] = 2 * func[i];
		deriv2[i] = 4 * func[i];
	}

	std::vector<double> deriv_fd = deriv_operator.mat_vec_prod(func);

	if (screen) {

		std::cout << std::scientific << std::setprecision(5);

//		std::cout << "Test derivative of exp function:" << std::endl << "Derivative order = " << deriv_order << std::endl;

//		std::cout << "dx    max-norm   l2-norm" << std::endl;

	}

	if (deriv_order == 1) {

		std::vector<double> diff = vector_diff(deriv1, deriv_fd);
		double max_norm_ = max_norm(diff);
		double dx = grid[1] - grid[0];
		double l2_norm_ = l2_norm(dx, diff);

		if (screen) {
//			print_test(grid, func, deriv1, deriv_fd);
			std::cout << std::setw(13) << dx << std::setw(13) << max_norm_ << std::setw(13) << l2_norm_ << std::endl;
		}

	}
	else if (deriv_order == 2) {

		std::vector<double> diff = vector_diff(deriv2, deriv_fd);
		double max_norm_ = max_norm(diff);
		double dx = grid[1] - grid[0];
		double l2_norm_ = l2_norm(dx, diff);

		if (screen) {
//			print_test(grid, func, deriv2, deriv_fd);
			std::cout << std::setw(13) << dx << std::setw(13) << max_norm_ << std::setw(13) << l2_norm_ << std::endl;
		}

	}

	return deriv_fd;

}


void test_cos(const std::vector<double>& grid, const int deriv_order, BandDiagonal& deriv_operator) {

	const int order = grid.size();

	std::vector<double> func(order, 0.0);
	std::vector<double> deriv1(order, 0.0);
	std::vector<double> deriv2(order, 0.0);

	for (int i = 0; i != order; ++i) {
		func[i] = cos(2 * grid[i]);
		deriv1[i] = -2 * sin(2 * grid[i]);
		deriv2[i] = -4 * func[i];
	}

	std::vector<double> deriv_fd = deriv_operator.mat_vec_prod(func);

	std::cout << "Test derivative of exp function:" << std::endl << "Derivative order = " << deriv_order << std::endl;

	if (deriv_order == 1) {
		print_test(grid, func, deriv1, deriv_fd);
	}
	else if (deriv_order == 2) {
		print_test(grid, func, deriv2, deriv_fd);
	}

}


void print_test(
	const std::vector<double>& grid,
	const std::vector<double>& func,
	const std::vector<double>& deriv,
	const std::vector<double>& deriv_fd) {

	std::cout << std::scientific << std::setprecision(5);

	std::cout << "       Grid           Function       Derivative     FD derivative  Abs diff" << std::endl;

	for (int i = 0; i != grid.size(); ++i) {
		std::cout
			<< std::setw(3) << i
			<< std::setw(15) << grid[i]
			<< std::setw(15) << func[i]
			<< std::setw(15) << deriv[i]
			<< std::setw(15) << deriv_fd[i]
			<< std::setw(15) << abs(deriv[i] - deriv_fd[i])
			<< std::endl;
	}

	std::cout << std::endl;

}


// Element-wise subtraction of vectors.
std::vector<double> vector_diff(
	const std::vector<double>& vec1,
	const std::vector<double>& vec2) {

	std::vector<double> diff(vec1.size(), 0.0);

	std::transform(vec1.begin(), vec1.end(), vec2.begin(), diff.begin(), std::minus<double>());

	return diff;

}


// Max norm.
double max_norm(std::vector<double> vec) {

	// Absolue value of each element.
	std::transform(vec.begin(), vec.end(), vec.begin(), [](double x) { return abs(x); });

	return *std::max_element(vec.begin(), vec.end());

}


// l2-norm (vector norm). TODO: Should this be dependent on dx? Like a "Riemann sum" expression?
double l2_norm(const double dx, std::vector<double> vec) {

	// Square of each element.
	std::transform(vec.begin(), vec.end(), vec.begin(), [dx](double x) { return dx * pow(x, 2); });

	return sqrt(std::accumulate(vec.begin(), vec.end(), (double)0));

}
