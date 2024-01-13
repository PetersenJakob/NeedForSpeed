#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "band_diagonal_matrix.h"


BandDiagonal::BandDiagonal(
	const int _order,
	const int _lower_bandwidth,
	const int _upper_bandwidth,
	const int _n_boundary_rows,
	const int _n_boundary_elements) {

	order_ = _order;
	lower_bandwidth_ = _lower_bandwidth;
	upper_bandwidth_ = _upper_bandwidth;
	bandwidth_ = 1 + lower_bandwidth_ + upper_bandwidth_;
	n_boundary_rows_ = _n_boundary_rows;
	n_boundary_elements_ = _n_boundary_elements;

	// First index: Diagonal index (sub -> main -> super).
	// Second index: Row index.
	std::vector<double> diagonal(order_, 0.0);
	std::vector<std::vector<double>> m(bandwidth_, diagonal);
	matrix = m;

	std::vector<double> row(n_boundary_elements_);
	std::vector<std::vector<double>> b(2 * n_boundary_rows_, row);
	boundary_rows = b;

}

// Adjust matrix rows at boundary using Gauss elimination.
void TriDiagonal::adjust_boundary(std::vector<double>& column) {

	if (n_boundary_rows_ != 1) {
		throw std::invalid_argument("Number of boundary rows should be 1.");
	}

	// TODO: what if n_boundary_elements_ != 3 and 4?

	if (n_boundary_elements_ == 2) {

		// Initialize "corner" elements not part of matrix.
		matrix[0][0] = 0.0;
		matrix[2][order_ - 1] = 0.0;

		for (int i = 0; i != n_boundary_elements_; ++i) {

			// Lower boundary row...
			matrix[i + 1][0] = boundary_rows[0][i];

			// Upper boundary row...
			matrix[i][order_ - 1] = boundary_rows[1][i];

		}
	}
	else if (n_boundary_elements_ == 3) {

		const double lower = boundary_rows[0][2] / matrix[2][1];
		const double upper = boundary_rows[1][0] / matrix[0][order_ - 2];

		for (int i = 0; i != n_boundary_elements_ - 1; ++i) {

			// Lower boundary row...
			matrix[i + 1][0] = boundary_rows[0][i] - lower * matrix[i][1];

			// Upper boundary row...
			matrix[i][order_ - 1] = boundary_rows[1][i + 1] - upper * matrix[i + 1][order_ - 2];

			// TODO: Adjust column vector on RHS of equal sign...
			column[0] -= lower * column[1];
			column[order_ - 1] -= upper * column[order_ - 2];

		}

		// Initialize "corner" elements not part of matrix.
		matrix[0][0] = 0.0;
		matrix[2][order_ - 1] = 0.0;

	}
	else if (n_boundary_elements_ == 4) {

		// Temporary boundary rows.
		std::vector<std::vector<double>> b_rows = boundary_rows;

		// First elimination.
		const double lower1 = b_rows[0][3] / matrix[2][2];
		const double upper1 = b_rows[1][0] / matrix[0][order_ - 3];

		for (int i = 0; i != 3; ++i) {

			// Lower boundary row...
			b_rows[0][i + 1] -= lower1 * matrix[i][2];

			// Upper boundary row...
			b_rows[1][i] -= upper1 * matrix[i][order_ - 3];

			// TODO: Adjust column vector on RHS of equal sign...

		}

		// Second elimination.
		const double lower2 = b_rows[0][2] / matrix[2][1];
		const double upper2 = b_rows[1][1] / matrix[0][order_ - 2];

		for (int i = 0; i != 2; ++i) {

			// Lower boundary row...
			matrix[i + 1][0] = b_rows[0][i] - lower2 * matrix[i][1];

			// Upper boundary row...
			matrix[i][order_ - 1] = b_rows[1][i + 2] - upper2 * matrix[i + 1][order_ - 2];

			// TODO: Adjust column vector on RHS of equal sign...

		}

		// Initialize "corner" elements not part of matrix.
		matrix[0][0] = 0.0;
		matrix[2][order_ - 1] = 0.0;

	}
}

std::vector<double> TriDiagonal::mat_vec_product(const std::vector<double>& column) {

	std::vector<double> result(order_, 0.0);

	for (int i = 0; i != n_boundary_rows_; ++i) {
		for (int j = 0; j != n_boundary_elements_; ++j) {
			result[i] += boundary_rows[i][j] * column[i + j];
		}
	}

	for (int i = n_boundary_rows_; i != order_ - n_boundary_rows_; ++i) {
		for (int j = 0; j != bandwidth_; ++j) {
			result[i] += matrix[j][i] * column[(i - n_boundary_rows_) + j];
		}
	}

	for (int i = n_boundary_rows_; i != 2 * n_boundary_rows_; ++i) {

		int row_nr = order_ - (2 * n_boundary_rows_ - i);

		for (int j = 0; j != n_boundary_elements_; ++j) {
			result[row_nr] += boundary_rows[i][j] * column[row_nr - (n_boundary_elements_ - 1) + j];
		}
	}

	return result;

}

// TODO: Why not const or reference allowed?
void print_matrix(BandDiagonal matrix) {

	std::cout << std::scientific << std::setprecision(5);

	std::cout
		<< "Band-diagonal matrix of order " << matrix.order() << ":" << std::endl
		<< "Lower bandwidth = " << matrix.lower_bandwidth() << std::endl
		<< "Upper bandwidth = " << matrix.upper_bandwidth() << std::endl;

	std::cout << "Boundary rows" << std::endl;
	for (int i = 0; i != matrix.n_boundary_rows(); ++i) {
		for (int j = 0; j != matrix.n_boundary_elements(); ++j) {
			std::cout << std::setw(14) << matrix.boundary_rows[i][j];
		}
		std::cout << std::endl;
	}

	std::cout << "Matrix" << std::endl;
	for (int i = 0; i != matrix.order(); ++i) {
		for (int j = 0; j != matrix.bandwidth(); ++j) {
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
