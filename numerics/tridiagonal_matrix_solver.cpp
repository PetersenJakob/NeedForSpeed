#include <vector>

#include "band_diagonal_matrix.h"
#include "tridiagonal_matrix_solver.h"

void tri_solver(TriDiagonal& matrix, std::vector<double>& column) {
	
//	matrix.adjust_boundary
	matrix.adjust_boundary(column);

	// Number of elements along main diagonal.
	const int n_elements = column.size();

	std::vector<double> tmp(n_elements, 0.0);

	// Forward sweep.
	double denominator = matrix.matrix[1][0];

	tmp[0] = matrix.matrix[2][0] / denominator;
	column[0] /= denominator;

	for (int i = 1; i != n_elements - 1; ++i) {

		denominator = matrix.matrix[1][i] - matrix.matrix[0][i] * tmp[i - 1];

		tmp[i] = matrix.matrix[2][i] / denominator;

		column[i] = (column[i] - matrix.matrix[0][i] * column[i - 1]) / denominator;

	}

	// Back substitution.
	for (int i = n_elements - 2; i != -1; --i) {
		column[i] -= tmp[i] * column[i + 1];
	}

}
