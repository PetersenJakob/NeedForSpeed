#pragma once

#include <functional>
#include <stdexcept>
#include <vector>

#include "matrix_equation_solver.h"


// Setting up finite difference representation of derivative operator on uniform grid.
template <class T>
T setup(
	const int order,
	const std::vector<double>& coef,
	const int n_boundary_rows,
	const int n_boundary_elements) {

	T matrix(order, n_boundary_rows, (n_boundary_rows - 1) + n_boundary_elements);

	for (int i = 0; i != coef.size(); ++i) {
		for (int j = n_boundary_rows; j != order - n_boundary_rows; ++j) {
			matrix.matrix[i][j] = coef[i];
		}
	}

	return matrix;

}


// Setting up finite difference representation of derivative operator on non-uniform grid.
template <class T>
T setup(
	const int order,
	const std::vector<double> grid,
	std::function<std::vector<double>(std::vector<double>)> coef,
	const int n_boundary_rows,
	const int n_boundary_elements) {

	T matrix(order, n_boundary_rows, (n_boundary_rows - 1) + n_boundary_elements);

	std::vector<double> coef_tmp;

	for (int i = n_boundary_rows; i != order - n_boundary_rows; ++i) {

		if (std::is_same<T, TriDiagonal>::value) {

			int idx_m1 = i - 1;
			int idx_p1 = i + 1;

			double dx_m1 = grid[i] - grid[idx_m1];
			double dx_p1 = grid[idx_p1] - grid[i];

			coef_tmp = coef({ 0.0, dx_m1, dx_p1, 0.0 });

		}
		else if (std::is_same<T, PentaDiagonal>::value) {

			int idx_m2 = i - 2;
			int idx_m1 = i - 1;
			int idx_p1 = i + 1;
			int idx_p2 = i + 2;

			double dx_m2 = grid[idx_m1] - grid[idx_m2];
			double dx_m1 = grid[i] - grid[idx_m1];
			double dx_p1 = grid[idx_p1] - grid[i];
			double dx_p2 = grid[idx_p2] - grid[idx_p1];

			coef_tmp = coef({ dx_m2, dx_m1, dx_p1, dx_p2 });

		}
		else {
			throw std::invalid_argument("Unknown matrix.");
		}

		for (int j = 0; j != matrix.n_diagonals(); ++j) {
			matrix.matrix[j][i] = coef_tmp[j];
		}

	}

	return matrix;

}


// Adjusting finite difference representation at boundary.
template <class T>
void boundary(const int row_index, const std::vector<double>& coef, T& matrix) {

	int index_tmp = row_index - matrix.n_boundary_rows();
	if (index_tmp < 0) {
		index_tmp = row_index;
	}

	for (int i = 0; i != coef.size(); ++i) {
		matrix.boundary_rows[row_index][index_tmp + i] = coef[i];
	}

}


// Evaulation of "derivative" operator wrt. first coordinate, 2-dimensional.
// solve_equation = true: Solve derivative * x = func
// solve_equation = false: Evaluate derivative * func
template <class T>
std::vector<double> action_2d(
	const int n_points_1,
	const int n_points_2,
	const int increment,
	const bool solve_equation,
	T& derivative,
	const std::vector<double>& func) {

	int factor_i = 1;
	int factor_j = 1;
	if (increment == 1) {
		// Function strip along x-dimension. TODO: (x, y)!
		factor_i = 1;
		factor_j = n_points_2;
	}
	else if (increment == 2) {
		// Function strip along y-dimension. TODO: (y, x)!
		factor_i = n_points_1;
		factor_j = 1;
	}
	else {
		throw std::invalid_argument("Unknown increment.");
	}

	std::vector<double> func_return(n_points_1 * n_points_2, 0.0);

	std::vector<double> func_strip(n_points_1, 0.0);

	int index = 0;

	for (int i = 0; i != n_points_2; ++i) {

		// Function strip along 1st dimension.
		for (int j = 0; j != n_points_1; ++j) {
			index = factor_i * i + factor_j * j;
			func_strip[j] = func[index];
		}

		// Evaluate expression.
		if (solve_equation) {
			solver::band(derivative, func_strip);
		}
		else {
			func_strip = derivative * func_strip;
		}

		// Save result.
		for (int j = 0; j != n_points_1; ++j) {
			index = factor_i * i + factor_j * j;
			func_return[index] = func_strip[j];
		}
	}

	return func_return;

}


// Evaulation of "derivative" operator wrt. first coordinate, 3-dimensional.
// solve_equation = true: Solve derivative * x = func
// solve_equation = false: Evaluate derivative * func
template <class T>
std::vector<double> action_3d(
	const int n_1,
	const int n_2,
	const int n_3,
	const int increment,
	const bool solve_equation,
	T& derivative,
	const std::vector<double>& func) {

	int factor_i = 1;
	int factor_j = 1;
	int factor_k = 1;
	if (increment == 1) {
		// Function strip along x-dimension. TODO: (x, y, z)!
		factor_i = n_3;
		factor_j = 1;
		factor_k = n_2 * n_3;
	}
	else if (increment == 2) {
		// Function strip along y-dimension. TODO: (y, x, z)!
		factor_i = n_1 * n_3;
		factor_j = 1;
		factor_k = n_3;
	}
	else if (increment == 3) {
		// Function strip along z-dimension. TODO: (z, x, y)!
		factor_i = n_1 * n_3;
		factor_j = n_1;
		factor_k = 1;
	}
	else {
		throw std::invalid_argument("Unknown increment.");
	}

	std::vector<double> func_result(n_1 * n_2 * n_3, 0.0);

	std::vector<double> func_strip(n_1, 0.0);

	int index = 0;

	for (int i = 0; i != n_2; ++i) {
		for (int j = 0; j != n_3; ++j) {

			// Function strip along 1st dimension.
			for (int k = 0; k != n_1; ++k) {
				index = factor_i * i + factor_j * j + factor_k * k;
				func_strip[k] = func[index];
			}

			// Evaluate expression.
			if (solve_equation) {
				solver::band(derivative, func_strip);
			}
			else {
				func_strip = derivative * func_strip;
			}

			for (int k = 0; k != n_1; ++k) {
				index = factor_i * i + factor_j * j + factor_k * k;
				func_result[index] = func_strip[k];
			}
		}
	}

	return func_result;

}
