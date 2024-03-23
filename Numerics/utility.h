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


// Evaulation of differential operator expression, 2-dimensional.
// Differential operator is wrt. first coordinate ("n_points_1").
// solve_equation
//	- true: differential * x = func
//  - false: x = differential * func
// Assume order of func to be (x, y).
template <class T>
std::vector<double> action_2d(
	const int n_points_1,
	const int n_points_2,
	const int filter,
	const bool solve_equation,
	T& derivative,
	const std::vector<double>& func) {

	int factor_i = 1;
	int factor_j = 1;

	if (filter == 1) {
		// Function strip along x-dimension. Order (x, y).
		factor_i = 1;
		factor_j = n_points_2;
	}
	else if (filter == 2) {
		// Function strip along y-dimension. Order (y, x).
		factor_i = n_points_1;
		factor_j = 1;
	}
	else {
		throw std::invalid_argument("Unknown filter.");
	}

	std::vector<double> func_strip(n_points_1, 0.0);

	const int n_points = n_points_1 * n_points_2;

	std::vector<double> func_return(n_points, 0.0);

	int index = 0;

	for (int i = 0; i != n_points_2; ++i) {

		// Function strip along 1st dimension.
		for (int j = 0; j != n_points_1; ++j) {
			index = factor_i * i + factor_j * j;
			func_strip[j] = func[index];
		}

		// Evaluate differential operator expression.
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


// Evaulation of differential operator expression, 2-dimensional.
// Differential operator is wrt. first coordinate ("n_points_1").
// solve_equation
//	- true: differential * x = func
//  - false: x = differential * func
// Assume order of func to be (x, y).
template <class T>
std::vector<double> action_2d(
	const int n_points_1,
	const int n_points_2,
	const int filter,
	const bool solve_equation,
	const std::vector<double>& adi_factors,
	const std::vector<double>& prefactors,
	std::vector<T>& derivatives,
	const std::vector<double>& func) {

	int factor_i = 1;
	int factor_j = 1;

	if (filter == 1) {
		// Function strip along x-dimension. Order (x, y).
		factor_i = 1;
		factor_j = n_points_2;
	}
	else if (filter == 2) {
		// Function strip along y-dimension. Order (y, x).
		factor_i = n_points_1;
		factor_j = 1;
	}
	else {
		throw std::invalid_argument("Unknown filter.");
	}

	std::vector<double> func_strip(n_points_1, 0.0);

	const int n_points = n_points_1 * n_points_2;

	std::vector<double> func_return(n_points, 0.0);

	int index = 0;

	T derivative = adi_factors[0] * derivatives[0]
		+ adi_factors[1] * (prefactors[0] * derivatives[1]
			+ prefactors[1] * derivatives[2]);

	for (int i = 0; i != n_points_2; ++i) {

		// Function strip along 1st dimension.
		for (int j = 0; j != n_points_1; ++j) {
			index = factor_i * i + factor_j * j;
			func_strip[j] = func[index];
		}

		// Evaluate differential operator expression.
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


// Evaulation of differential operator expression, 2-dimensional.
// Differential operator is wrt. first coordinate ("n_points_1").
// solve_equation
//	- true: differential * x = func
//  - false: x = differential * func
// Assume order of func to be (x, y).
template <class T>
std::vector<double> action_2d(
	const int n_points_1,
	const int n_points_2,
	const int filter,
	const bool solve_equation,
	const std::vector<double> adi_factors,
	const std::vector<std::vector<double>>& prefactors,
	std::vector<T>& derivatives,
	const std::vector<double>& func) {

	int factor_i = 1;
	int factor_j = 1;

	if (filter == 1) {
		// Function strip along x-dimension. Order (x, y).
		factor_i = 1;
		factor_j = n_points_2;
	}
	else if (filter == 2) {
		// Function strip along y-dimension. Order (y, x).
		factor_i = n_points_1;
		factor_j = 1;
	}
	else {
		throw std::invalid_argument("Unknown filter.");
	}

	std::vector<double> func_strip(n_points_1, 0.0);

	const int n_points = n_points_1 * n_points_2;

	std::vector<double> func_return(n_points, 0.0);

	int index = 0;



	double factor1 = prefactors[0][0];
	double factor2 = prefactors[1][0];

	double factor3 = prefactors[2][0];

	int n_start = 0;
	int n_final = 0;
	int n_index = 0;

	if (filter == 1) {
		n_start = 1;
		n_final = 1 + n_points_1;
		n_index = 1 + n_points_1;
	}
	else {
		n_start = 1 + n_points_2;
		n_final = 1 + n_points_2 + n_points_1;
		n_index = 1;
	}

	std::vector<double> vec1 =
		std::vector<double>(prefactors[0].begin() + n_start, prefactors[0].begin() + n_final);
	std::vector<double> vec2 =
		std::vector<double>(prefactors[1].begin() + n_start, prefactors[1].begin() + n_final);

	std::vector<double> vec3 =
		std::vector<double>(prefactors[2].begin() + n_start, prefactors[2].begin() + n_final);

	for (int j = 0; j != n_points_1; ++j) {
		vec1[j] *= factor1;
		vec2[j] *= factor2;

		vec3[j] *= factor3;

	}

	std::vector<double> vec1_tmp = vec1;
	std::vector<double> vec2_tmp = vec2;

	std::vector<double> vec3_tmp = vec3;


	T derivative = derivatives[0];
	T deriv_0 = adi_factors[0] * derivatives[0];
	T deriv_1 = adi_factors[1] * derivatives[1];
	T deriv_2 = adi_factors[1] * derivatives[2];

	T inhomo = adi_factors[1] * derivatives[0];


	for (int i = 0; i != n_points_2; ++i) {

		// Function strip along 1st dimension.
		for (int j = 0; j != n_points_1; ++j) {
			index = factor_i * i + factor_j * j;
			func_strip[j] = func[index];
		}

		// Update derivative operator.
		// identity + c1 * f1(x) * g1(x) * d1dx1 + c2 * f2(y) * g2(y) * d2dx2 + c3 * f3(x) * g3(y) * identity.
		for (int j = 0; j != n_points_1; ++j) {
			index = n_index + i;
			vec1_tmp[j] = vec1[j] * prefactors[0][index];
			vec2_tmp[j] = vec2[j] * prefactors[1][index];

			vec3_tmp[j] = vec3[j] * prefactors[2][index];

		}

		derivative = deriv_0;
		derivative += deriv_1.pre_vector(vec1_tmp);
 		derivative += deriv_2.pre_vector(vec2_tmp);

		derivative += inhomo.pre_vector(vec3_tmp);

		// Evaluate differential operator expression.
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


// Evaulation of differential operator expression, 3-dimensional.
// Differential operator is wrt. first coordinate ("n_points_1").
// solve_equation
//	- true: differential * x = func
//  - false: x = differential * func
// Assume order of func to be (x, y, z).
template <class T>
std::vector<double> action_3d(
	const int n_points_1,
	const int n_points_2,
	const int n_points_3,
	const int filter,
	const bool solve_equation,
	T& derivative,
	const std::vector<double>& func) {

	int factor_i = 1;
	int factor_j = 1;
	int factor_k = 1;

	if (filter == 1) {
		// Function strip along x-dimension. Order (x, y, z).
		factor_i = n_points_3;
		factor_j = 1;
		factor_k = n_points_2 * n_points_3;
	}
	else if (filter == 2) {
		// Function strip along y-dimension. Order (y, x, z).
		factor_i = n_points_1 * n_points_3;
		factor_j = 1;
		factor_k = n_points_3;
	}
	else if (filter == 3) {
		// Function strip along z-dimension. Order (z, x, y).
		factor_i = n_points_1 * n_points_3;
		factor_j = n_points_1;
		factor_k = 1;
	}
	else {
		throw std::invalid_argument("Unknown filter.");
	}

	std::vector<double> func_strip(n_points_1, 0.0);

	const int n_points = n_points_1 * n_points_2 * n_points_3;

	std::vector<double> func_result(n_points, 0.0);

	int index = 0;

	for (int i = 0; i != n_points_2; ++i) {

		for (int j = 0; j != n_points_3; ++j) {

			// Function strip along 1st dimension.
			for (int k = 0; k != n_points_1; ++k) {
				index = factor_i * i + factor_j * j + factor_k * k;
				func_strip[k] = func[index];
			}

			// Evaluate differential operator expression.
			if (solve_equation) {
				solver::band(derivative, func_strip);
			}
			else {
				func_strip = derivative * func_strip;
			}

			// Save result.
			for (int k = 0; k != n_points_1; ++k) {
				index = factor_i * i + factor_j * j + factor_k * k;
				func_result[index] = func_strip[k];
			}

		}

	}

	return func_result;

}


// Evaulation of differential operator expression, 4-dimensional.
// Differential operator is wrt. first coordinate ("n_points_1").
// solve_equation
//	- true: differential * x = func
//  - false: x = differential * func
// Assume order of func to be (x1, x2, x3, x4).
template <class T>
std::vector<double> action_4d(
	const int n_points_1,
	const int n_points_2,
	const int n_points_3,
	const int n_points_4,
	const int filter,
	const bool solve_equation,
	T& derivative,
	const std::vector<double>& func) {

	int factor_i = 1;
	int factor_j = 1;
	int factor_k = 1;
	int factor_l = 1;


	if (filter == 1) {
		// Function strip along x1-dimension. Order (x1, x2, x3, x4).
		factor_i = n_points_3 * n_points_4;
		factor_j = n_points_4;
		factor_k = 1;
		factor_l = n_points_2 * n_points_3 * n_points_4;
	}
	else if (filter == 2) {
		// Function strip along x2-dimension. Order (x2, x1, x3, x4).
		factor_i = n_points_1 * n_points_3 * n_points_4;
		factor_j = n_points_4;
		factor_k = 1;
		factor_l = n_points_3 * n_points_4;
	}
	else if (filter == 3) {
		// Function strip along x3-dimension. Order (x3, x1, x2, x4).
		factor_i = n_points_1 * n_points_3 * n_points_4;
		factor_j = n_points_1 * n_points_4;
		factor_k = 1;
		factor_l = n_points_4;
	}
	else if (filter == 4) {
		// Function strip along x4-dimension. Order (x4, x1, x2, x3).
		factor_i = n_points_1 * n_points_3 * n_points_4;
		factor_j = n_points_1 * n_points_4;
		factor_k = n_points_1;
		factor_l = 1;
	}
	else {
		throw std::invalid_argument("Unknown filter.");
	}

	std::vector<double> func_strip(n_points_1, 0.0);

	const int n_points = n_points_1 * n_points_2 * n_points_3 * n_points_4;

	std::vector<double> func_result(n_points, 0.0);

	int index = 0;

	for (int i = 0; i != n_points_2; ++i) {

		for (int j = 0; j != n_points_3; ++j) {

			for (int k = 0; k != n_points_4; ++k) {

				// Function strip along 1st dimension.
				for (int l = 0; l != n_points_1; ++l) {
					index = 
						factor_i * i + factor_j * j + factor_k * k + factor_l * l;
					func_strip[l] = func[index];
				}

				// Evaluate differential operator expression.
				if (solve_equation) {
					solver::band(derivative, func_strip);
				}
				else {
					func_strip = derivative * func_strip;
				}

				// Save result.
				for (int l = 0; l != n_points_1; ++l) {
					index = 
						factor_i * i + factor_j * j + factor_k * k + factor_l * l;
					func_result[index] = func_strip[l];
				}

			}

		}

	}

	return func_result;

}
