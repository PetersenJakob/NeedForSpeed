#include <iostream>
#include <vector>

#include "band_diagonal_matrix.h"
#include "tridiagonal_matrix_solver.h"

void tri_solver(TriDiagonal& matrix, std::vector<double>& column) {
	
#if false
	print_matrix(matrix);
	std::cout << "Column vector before: " << std::endl;
	for (int i = 0; i != column.size(); ++i) {
		std::cout << column[i] << std::endl;
	}
	std::cout << std::endl;
#endif
//	matrix.adjust_boundary
	matrix.adjust_boundary(column);
#if false
	print_matrix(matrix);
	std::cout << "Column vector after: " << std::endl;
	for (int i = 0; i != column.size(); ++i) {
		std::cout << column[i] << std::endl;
	}
	std::cout << std::endl;
#endif
	// Number of elements along main diagonal.
	const int n_elements = column.size();

	std::vector<double> tmp(n_elements, 0.0);

	// Forward sweep.
	double denominator = matrix.matrix[1][0];

	tmp[0] = matrix.matrix[2][0] / denominator;
	column[0] /= denominator;

	int idx_tmp;
	for (int i = 1; i != n_elements; ++i) {

		idx_tmp = i - 1;

		denominator = matrix.matrix[1][i] - matrix.matrix[0][i] * tmp[idx_tmp];

		tmp[i] = matrix.matrix[2][i] / denominator;

		column[i] = (column[i] - matrix.matrix[0][i] * column[idx_tmp]) / denominator;

	}

	// Back substitution.
	for (int i = n_elements - 2; i != -1; --i) {

		idx_tmp = i + 1;

		column[i] -= tmp[i] * column[idx_tmp];

	}

}
