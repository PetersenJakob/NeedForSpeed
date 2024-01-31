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

	std::vector<double> vec_tmp(matrix.order(), 0.0);

	tridiagonal_matrix_solver(
		matrix.matrix[0],
		matrix.matrix[1],
		matrix.matrix[2],
		column,
		vec_tmp);

}


void tridiagonal_matrix_solver(
	std::vector<double>& sub,
	std::vector<double>& main,
	std::vector<double>& super,
	std::vector<double>& column,
	std::vector<double>& vec_tmp) {

	// Number of elements along main diagonal.
	const int n_elements = main.size();

	// Temporary index.
	int idx_tmp = 0;

	// *****************************************************
	// Forward sweep:
	// Remove elements of sub-diagonal by Gauss elimination.
	// *****************************************************

	// Denominator for normalization of 1st element of main diagonal.
	double denominator = main[0];

	// 1st element of super-diagonal after normalization.
	vec_tmp[0] = super[0] / denominator;

	// 1st element of column vector after normalization.
	column[0] /= denominator;

	for (int i = 1; i != n_elements; ++i) {

		idx_tmp = i - 1;

		// Denominator for normalization of (i + 1)'th element of main diagonal.
		denominator = main[i] - sub[i] * vec_tmp[idx_tmp];

		// (i + 1)'th element of super-diagonal after Gauss elimination.
		vec_tmp[i] = super[i] / denominator;

		// (i + 1)'th element of column vector after Gauss elimination.
		column[i] = (column[i] - sub[i] * column[idx_tmp]) / denominator;

	}

	// ******************************************************
	// Back substitution:
	// Remove element of super-diagonal by Gauss elimination.
	// ******************************************************

	for (int i = n_elements - 2; i != -1; --i) {

		idx_tmp = i + 1;

		// (i + 1)'th element of vector matrix after Gauss elimination.
		column[i] -= vec_tmp[i] * column[idx_tmp];

	}

}


void pentadiagonal_matrix_solver(
	std::vector<double>& sub_2,
	std::vector<double>& sub_1,
	std::vector<double>& main,
	std::vector<double>& super_1,
	std::vector<double>& super_2,
	std::vector<double>& column,
	std::vector<double>& sub_tmp,
	std::vector<double>& main_tmp,
	std::vector<double>& super_tmp,
	std::vector<double>& vec_tmp) {

	// Number of elements along main diagonal.
	const int n_elements = main.size();

	// Temporary index.
	int idx_tmp = 0;

	// *****************************************************
	// Forward sweep:
	// Remove elements of sub-diagonal by Gauss elimination.
	// *****************************************************

	// Denominator for normalization of 1st element of main diagonal.
	double denominator = main[0];

	// 1st element of super-diagonal after normalization.
	vec_tmp[0] = super_1[0] / denominator;

	// 1st element of column vector after normalization.
	column[0] /= denominator;

	for (int i = 1; i != n_elements; ++i) {

		idx_tmp = i - 1;

		// Denominator for normalization of (i + 1)'th element of main diagonal.
		denominator = main[i] - sub_1[i] * vec_tmp[idx_tmp];

		// (i + 1)'th element of super-diagonal after Gauss elimination.
		vec_tmp[i] = super_1[i] / denominator;

		// (i + 1)'th element of column vector after Gauss elimination.
		column[i] = (column[i] - sub_1[i] * column[idx_tmp]) / denominator;

	}

	// ******************************************************
	// Back substitution:
	// Remove element of super-diagonal by Gauss elimination.
	// ******************************************************

	for (int i = n_elements - 2; i != -1; --i) {

		idx_tmp = i + 1;

		// (i + 1)'th element of vector matrix after Gauss elimination.
		column[i] -= vec_tmp[i] * column[idx_tmp];

	}

	// ...
	tridiagonal_matrix_solver(sub_tmp, main_tmp, super_tmp, column, vec_tmp);

}
