#define _USE_MATH_DEFINES
#include <cmath>
#include <iomanip>
#include <iostream>

#include "band_diagonal_matrix.h"
#include "finite_difference.h"
#include "grid_generator.h"
#include "tridiagonal_matrix_solver.h"


int main() {

	const double x_min = -0.2; // -M_PI / 2.0;
	const double x_max = 0.2; // M_PI / 2.0;
	const int order = 21;
	const double dx = (x_max - x_min) / (order - 1);

	const double dt = 0.1;

	std::vector<double> column_tmp(order, 0.0);

	TriDiagonal tri(order);
//	print_matrix(tri);

	PentaDiagonal penta(order);
//	print_matrix(penta);


	PentaDiagonal d1dx1_p = d1dx1::c4b4(order, dx);
	print_matrix(d1dx1_p);

	d1dx1_p.adjust_boundary(column_tmp);
	print_matrix(d1dx1_p);


	TriDiagonal d1dx1_ = d1dx1::c2b2(order, dx);
//	print_matrix(d1dx1_);

	d1dx1_.adjust_boundary(column_tmp);
//	print_matrix(d1dx1_);


	std::vector<double> grid = grid_equidistant(x_min, x_max, order);
	std::cout << std::scientific << std::setprecision(5);

	std::vector<double> func;
	for (int i = 0; i != order; ++i)
		func.push_back(exp(2 * grid[i]));

	std::vector<double> deriv1 = d1dx1_.mat_vec_prod(func);

	for (int i = 0; i != order; ++i) {

		std::cout 
			<< std::setw(3) << i 
			<< std::setw(14) << grid[i] 
			<< std::setw(14) << func[i]
			<< std::setw(14) << 2 * exp(2 * grid[i])
			<< std::setw(14) << deriv1[i]
			<< std::endl;

	}

	TriDiagonal d2dx2_ = d2dx2::c2b2(order, dx);
	print_matrix(d2dx2_);

	d2dx2_.adjust_boundary(column_tmp);
	print_matrix(d2dx2_);

	std::vector<double> deriv2 = d2dx2_.mat_vec_prod(func);

	for (int i = 0; i != order; ++i) {

		std::cout
			<< std::setw(3) << i
			<< std::setw(14) << grid[i]
			<< std::setw(14) << func[i]
			<< std::setw(14) << 4 * exp(2 * grid[i])
			<< std::setw(14) << deriv2[i]
			<< std::endl;

	}
	std::cout << std::endl;



	std::vector<double> function(order, 0.0);
	std::vector<double> antideriv(order, 0.0);
	std::vector<double> antiantideriv(order, 0.0);
	std::vector<double> deriv(order, 0.0);
	std::vector<double> derivderiv(order, 0.0);

	for (int i = 0; i != order; ++i) {
		function[i] = cos(grid[i]);
		antideriv[i] = sin(grid[i]);
		antiantideriv[i] = -function[i];
		deriv[i] = -sin(grid[i]);
		derivderiv[i] = -function[i];
	}

	std::vector<double> deriv_fd = d1dx1_.mat_vec_prod(function);
	std::vector<double> derivderiv_fd = d2dx2_.mat_vec_prod(function);

	std::vector<double> solution_fd = function;

	for (int i = 0; i != d2dx2_.order(); ++i) {

		for (int j = 0; j != d2dx2_.n_diagonals(); ++j) {
			d2dx2_.matrix[j][i] *= -dt;
		}

		d2dx2_.matrix[d2dx2_.bandwidth()][i] += 1.0;

	}

	tri_solver(d2dx2_, solution_fd);

	for (int i = 0; i != order; ++i) {

		std::cout
			<< std::setw(3) << i
			<< std::setw(14) << grid[i]
			<< std::setw(14) << function[i]
			<< std::setw(14) << deriv[i]
			<< std::setw(14) << deriv_fd[i]
			<< std::setw(14) << exp(-dt) * function[i]
			<< std::setw(14) << solution_fd[i]
			<< std::endl;

	}

	return 0;

}