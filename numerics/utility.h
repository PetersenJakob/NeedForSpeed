#pragma once

#include <functional>
#include <stdexcept>
#include <vector>


// TODO: Why should these template functions be defined in the header file?

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
