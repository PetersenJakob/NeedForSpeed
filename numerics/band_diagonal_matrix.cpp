#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "band_diagonal_matrix.h"


BandDiagonal::BandDiagonal(
	const int _order,
	const int _bandwidth,
	const int _n_boundary_rows,
	const int _n_boundary_elements) {

	order_ = _order;
	bandwidth_ = _bandwidth;
	n_diagonals_ = 1 + 2 * _bandwidth;
	n_boundary_rows_ = _n_boundary_rows;
	n_boundary_elements_ = _n_boundary_elements;

	// Band-diagonal matrix in compact form.
	std::vector<double> diagonal(order_, 0.0);
	std::vector<std::vector<double>> m(n_diagonals_, diagonal);
	matrix = m;

	// Boundary rows of band-diagonal matrix.
	std::vector<double> row(n_boundary_elements_);
	std::vector<std::vector<double>> b(2 * n_boundary_rows_, row);
	boundary_rows = b;

}


std::vector<double> BandDiagonal::mat_vec_prod(const std::vector<double>& column) {

	int row_idx;
	int col_idx;
	std::vector<double> result(order_, 0.0);

	// Lower boundary rows.
	for (int i = 0; i != n_boundary_rows_; ++i) {
		for (int j = 0; j != n_boundary_elements_; ++j) {
			col_idx = i + j;
			result[i] += boundary_rows[i][j] * column[col_idx];
		}
	}

	// Interior rows.
	for (int i = n_boundary_rows_; i != order_ - n_boundary_rows_; ++i) {
		for (int j = 0; j != n_diagonals_; ++j) {
			col_idx = (i - n_boundary_rows_) + j;
			result[i] += matrix[j][i] * column[col_idx];
		}
	}

	// Upper boundary rows.
	for (int i = n_boundary_rows_; i != 2 * n_boundary_rows_; ++i) {
		row_idx = (order_ - 2 * n_boundary_rows_) + i;
		for (int j = 0; j != n_boundary_elements_; ++j) {
			col_idx = (row_idx - n_boundary_elements_) + j + 1;
			result[row_idx] += boundary_rows[i][j] * column[col_idx];
		}
	}

	return result;

}


// Remove boundary row element by Gauss elimination.
void BandDiagonal::gauss_elimination(
	const int boundary_row_idx,
	const int boundary_element_idx,
	const int matrix_row_idx,
	std::vector<double>& column) {

	// Lower boundary row index.
	const int br_lower_idx = boundary_row_idx;
	// Lower boundary row element index. TODO: Remember to include zero at beginning if more than one boundary row!
	const int be_lower_idx = boundary_element_idx;
	// Lower matrix row index.
	const int mr_lower_idx = matrix_row_idx;
	// Lower matrix row element index. TODO: Assume that the last element is used to remove boundary element!
	const int me_lower_idx = n_diagonals_ - 1;

	// Upper boundary row index.
	const int br_upper_idx = (2 * n_boundary_rows_ - 1) - br_lower_idx;
	// Upper boundary row element index.
	const int be_upper_idx = (n_boundary_elements_ - 1) - be_lower_idx;
	// Upper matrix row index.
	const int mr_upper_idx = (order_ - 1) - mr_lower_idx;
	// Upper matrix row element index.
	const int me_upper_idx = (n_diagonals_ - 1) - me_lower_idx;

	// Factor at lower boundary.
	const double lower = boundary_rows_tmp[br_lower_idx][be_lower_idx] / matrix[me_lower_idx][mr_lower_idx];
	// Factor at upper boundary.
	const double upper = boundary_rows_tmp[br_upper_idx][be_upper_idx] / matrix[me_upper_idx][mr_upper_idx];

	int me_lower_idx_tmp = 0;
	int be_lower_idx_tmp = 0;
	int me_upper_idx_tmp = 0;
	int be_upper_idx_tmp = 0;

	for (int i = 0; i != n_diagonals_; ++i) {

		// Adjust lower boundary rows.
		me_lower_idx_tmp = me_lower_idx - i;
		be_lower_idx_tmp = be_lower_idx - i;
		boundary_rows_tmp[br_lower_idx][be_lower_idx_tmp] -= lower * matrix[me_lower_idx_tmp][mr_lower_idx];

		// Adjust upper boundary rows.
		me_upper_idx_tmp = me_upper_idx + i;
		be_upper_idx_tmp = be_upper_idx + i;
		boundary_rows_tmp[br_upper_idx][be_upper_idx_tmp] -= lower * matrix[me_upper_idx_tmp][mr_upper_idx];

		// Adjust RHS column vector.
		column[br_lower_idx] -= lower * column[mr_lower_idx];
		column[br_upper_idx] -= lower * column[mr_upper_idx];

	}

}


// Overwrite boundary row of matrix.
void BandDiagonal::overwrite_bounary_row(const int boundary_row_idx) {

	// Lower boundary row index.
	const int br_lower_idx = boundary_row_idx;
	// Lower matrix row index.
	const int mr_lower_idx = boundary_row_idx;

	// Upper boundary row index.
	const int br_upper_idx = (2 * n_boundary_rows_ - 1) - br_lower_idx;
	// Upper matrix row index.
	const int mr_upper_idx = (order_ - 1) - mr_lower_idx;

	for (int i = 0; i != boundary_row_idx + (bandwidth_ + 1); ++i) {

		int me_lower_idx = (bandwidth_ - boundary_row_idx) + i;
		int me_upper_idx = i;

		int be_lower_idx = i;
		int be_upper_idx = n_boundary_elements_ - (boundary_row_idx + (bandwidth_ + 1)) + i ;

		matrix[me_lower_idx][mr_lower_idx] = boundary_rows_tmp[br_lower_idx][be_lower_idx];

		matrix[me_upper_idx][mr_upper_idx] = boundary_rows_tmp[br_upper_idx][be_upper_idx];

	}

}


// Adjust matrix rows at boundary using Gauss elimination.
void TriDiagonal::adjust_boundary(std::vector<double>& column) {

	if (n_boundary_rows_ != 1) {
		throw std::invalid_argument("Number of boundary rows should be 1.");
	}

	// TODO: what if n_boundary_elements_ != 3 and 4?

	if (n_boundary_elements_ == 2) {

		boundary_rows_tmp = boundary_rows;
		overwrite_bounary_row(0);

	}
	else if (n_boundary_elements_ == 3) {

		boundary_rows_tmp = boundary_rows;
		gauss_elimination(0, 2, 1, column);
		overwrite_bounary_row(0);

	}
	else if (n_boundary_elements_ == 4) {

		boundary_rows_tmp = boundary_rows;
		gauss_elimination(0, 3, 2, column);
		gauss_elimination(0, 2, 1, column);
		overwrite_bounary_row(0);

	}
	else if (n_boundary_elements_ == 5) {

		boundary_rows_tmp = boundary_rows;
		gauss_elimination(0, 4, 3, column);
		gauss_elimination(0, 3, 2, column);
		gauss_elimination(0, 2, 1, column);
		overwrite_bounary_row(0);

	}

	// Initialize "corner" elements which are not part of matrix.
	matrix[0][0] = 0.0;
	matrix[2][order_ - 1] = 0.0;

}

// TODO: Why not const or reference allowed?
void print_matrix(BandDiagonal matrix) {

	std::cout << std::scientific << std::setprecision(5);

	std::cout << std::endl
		<< "Band-diagonal matrix:" << std::endl
		<< "Order = " << matrix.order()
		<< ", bandwidth = " << matrix.bandwidth()
		<< ", n_diagonals = " << matrix.n_diagonals()
		<< ", n_boundary_rows = " << matrix.n_boundary_rows()
		<< ", n_boundary_elements = " << matrix.n_boundary_elements() << std::endl;

	std::cout << "Boundary rows" << std::endl;
	for (int i = 0; i != matrix.n_boundary_rows(); ++i) {
		for (int j = 0; j != matrix.n_boundary_elements(); ++j) {
			std::cout << std::setw(14) << matrix.boundary_rows[i][j];
		}
		std::cout << std::endl;
	}

	std::cout << "Matrix" << std::endl;
	for (int i = 0; i != matrix.order(); ++i) {
		for (int j = 0; j != matrix.n_diagonals(); ++j) {
			std::cout << std::setw(14) << matrix.matrix[j][i];
		}
		std::cout << std::endl;
	}

	std::cout << "Boundary rows" << std::endl;
	for (int i = matrix.n_boundary_rows(); i != 2 * matrix.n_boundary_rows(); ++i) {
		for (int j = 0; j != matrix.n_boundary_elements(); ++j) {
			std::cout << std::setw(14) << matrix.boundary_rows[i][j];
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

}
