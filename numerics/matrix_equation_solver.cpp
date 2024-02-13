#include <iostream>
#include <vector>

#include "band_diagonal_matrix.h"
#include "matrix_equation_solver.h"


void tri_solver(TriDiagonal& matrix, std::vector<double>& column) {

	matrix.adjust_boundary(column);

	std::vector<double> vec_tmp(matrix.order(), 0.0);

	tridiagonal_matrix_solver(
		matrix.matrix[0],
		matrix.matrix[1],
		matrix.matrix[2],
		column,
		vec_tmp);

}


void penta_solver(PentaDiagonal& matrix, std::vector<double>& column) {

	matrix.adjust_boundary(column);

	std::vector<double> sub_tmp(matrix.order(), 0.0);
	std::vector<double> main_tmp(matrix.order(), 0.0);
	std::vector<double> super_tmp(matrix.order(), 0.0);
	std::vector<double> vec_tmp(matrix.order(), 0.0);

	pentadiagonal_matrix_solver(
		matrix.matrix[0],
		matrix.matrix[1],
		matrix.matrix[2],
		matrix.matrix[3],
		matrix.matrix[4],
		column,
		sub_tmp,
		main_tmp,
		super_tmp,
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

	// *******************************************************
	// Backward sweep (back substitution):
	// Remove elements of super-diagonal by Gauss elimination.
	// *******************************************************

	for (int i = n_elements - 2; i != -1; --i) {

		idx_tmp = i + 1;

		// (i + 1)'th element of column vector after Gauss elimination.
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

	// *********************************************************
	// Forward sweep:
	// Remove elements of 2nd sub-diagonal by Gauss elimination.
	// *********************************************************

	// Denominator for normalization of 1st element of main diagonal.
	double denominator = main[0];


	// 1st element of main diagonal after normalization.
	main_tmp[0] = 1.0;

	// 1st element of super-diagonals after normalization.
	super_tmp[0] = super_1[0] / denominator;
	vec_tmp[0] = super_2[0] / denominator;

	// 1st element of column vector after normalization.
	column[0] /= denominator;


	// Add 1st row to 2nd row. !!! Avoid zero as 2nd element of 1st sub-diagonal !!!
	sub_tmp[1] = sub_1[1] + main_tmp[0];
	main_tmp[1] = main[1] + super_tmp[0];
	super_tmp[1] = super_1[1] + vec_tmp[0];
	vec_tmp[1] = super_2[1];
	column[1] += column[0];

#if false
	sub_tmp[1] = sub_1[1];
	main_tmp[1] = main[1];
	super_tmp[1] = super_1[1];
	vec_tmp[1] = super_2[1];
#endif

	// Denominator for normalization of 2nd element of 1st sub-diagonal.
	denominator = sub_tmp[1];

	// 2nd element of 1st sub-diagonal after normalization.
	sub_tmp[1] = 1.0;

	// 2nd element of main diagonal after normalization.
	main_tmp[1] /= denominator;

	// 2nd element of super-diagonals after normalization.
	super_tmp[1] /= denominator;
	vec_tmp[1] /= denominator;

	// 2nd element of column vector after normalization.
	column[1] /= denominator;


	// TODO: The range depends on choice of boundary rows!
	for (int i = 2; i != n_elements - 2; ++i) {

		idx_tmp = i - 1;

		// Denominator for normalization of (i + 1)'th element of 1st sub-diagonal.
		denominator = sub_1[i] - sub_2[i] * main_tmp[idx_tmp];

		// (i + 1)'th element of 1st sub-diagonal after Gauss elimination.
		sub_tmp[i] = 1.0;

		// (i + 1)'th element of main diagonal after Gauss elimination.
		main_tmp[i] = (main[i] - sub_2[i] * super_tmp[idx_tmp]) / denominator;

		// (i + 1)'th element of super-diagonals after Gauss elimination.
		super_tmp[i] = (super_1[i] - sub_2[i] * vec_tmp[idx_tmp]) / denominator;
		vec_tmp[i] = super_2[i] / denominator;

		// (i + 1)'th element of column vector after Gauss elimination.
		column[i] = (column[i] - sub_2[i] * column[idx_tmp]) / denominator;

	}

	// ***********************************************************
	// Backward sweep:
	// Remove elements of 2nd super-diagonal by Gauss elimination.
	// ***********************************************************


	// Index of last and 2nd last element.
	const int idx_last = n_elements - 1;
	const int idx_2nd_last = n_elements - 2;

#if false
	// Denominator for normalization of last element of main diagonal.
	denominator = main_tmp[idx_last];

	// Last element of main diagonal after normalization.
	main_tmp[idx_last] = 1.0;

	// Last element of sub-diagonal after normalization.
	sub_tmp[idx_last] /= denominator;

	// Last element of column vector after normalization.
	column[idx_last] /= denominator;


	// Add last row to 2nd last row. !!! Avoid zero as 2nd element of 1st sub-diagonal !!!
	main_tmp[idx_2nd_last] += sub_tmp[idx_last];
	super_tmp[idx_2nd_last] += main_tmp[idx_last];
	column[idx_2nd_last] += column[idx_last];


	// Denominator for normalization of 2nd last element of 1st super-diagonal.
	denominator = super_tmp[idx_2nd_last];

	// 2nd last element of 1st super-diagonal after normalization.
	super_tmp[idx_2nd_last] = 1.0;

	// 2nd last element of main diagonal after normalization.
	main_tmp[idx_2nd_last] /= denominator;

	// 2nd last element of 1st sub-diagonal after normalization.
	sub_tmp[idx_2nd_last] /= denominator;

	// 2nd last element of column vector after normalization.
	column[idx_2nd_last] /= denominator;
#endif

	main_tmp[idx_last] = main[idx_last];

	sub_tmp[idx_2nd_last] = sub_1[idx_2nd_last];
	main_tmp[idx_2nd_last] = main[idx_2nd_last];
	super_tmp[idx_2nd_last] = super_1[idx_2nd_last];

	// Denominator for normalization of 2nd last element of 1st super-diagonal.
	denominator = super_tmp[idx_2nd_last];

	// 2nd last element of 1st super-diagonal after normalization.
	super_tmp[idx_2nd_last] = 1.0;

	// 2nd last element of main diagonal after normalization.
	main_tmp[idx_2nd_last] /= denominator;

	// 2nd last element of 1st sub-diagonal after normalization.
	sub_tmp[idx_2nd_last] /= denominator;

	// 2nd last element of column vector after normalization.
	column[idx_2nd_last] /= denominator;

//	for (int i = n_elements - 3; i != -1; --i) {
	for (int i = n_elements - 3; i != 1; --i) {

		idx_tmp = i + 1;

		// Denominator for normalization of (i + 1)'th element of 1st super-diagonal.
		denominator = super_tmp[i] - vec_tmp[i] * main_tmp[idx_tmp];

		// (i + 1)'th element of 1st super-diagonal after Gauss elimination.
		super_tmp[i] = 1.0;

		// (i + 1)'th element of main diagonal after Gauss elimination.
		main_tmp[i] = (main_tmp[i] - vec_tmp[i] * sub_tmp[idx_tmp]) / denominator;

		// (i + 1)'th element of 1st sub-diagonal after Gauss elimination.
		sub_tmp[i] /= denominator;

		// (i + 1)'th element of column vector after Gauss elimination.
		column[i] = (column[i] - vec_tmp[i] * column[idx_tmp]) / denominator;

	}

	// ...
	tridiagonal_matrix_solver(sub_tmp, main_tmp, super_tmp, column, vec_tmp);

}
