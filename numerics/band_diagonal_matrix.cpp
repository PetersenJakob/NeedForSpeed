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


BandDiagonal::BandDiagonal(const BandDiagonal& mat) {

	order_ = mat.order_;
	bandwidth_ = mat.bandwidth_;
	n_diagonals_ = mat.n_diagonals_;
	n_boundary_rows_ = mat.n_boundary_rows_;
	n_boundary_elements_ = mat.n_boundary_elements_;
	matrix = mat.matrix;
	boundary_rows = mat.boundary_rows;
	boundary_rows_tmp = mat.boundary_rows_tmp;

}


bool BandDiagonal::operator==(const BandDiagonal& m)
{
	const double eps = 1.0e-8;

	if (order_ == m.order_ &&
		bandwidth_ == m.bandwidth_ &&
		n_boundary_rows_ == m.n_boundary_rows_ &&
		n_boundary_elements_ == m.n_boundary_elements_) {

		for (int i = n_boundary_rows_; i != order_ - n_boundary_rows_; ++i) {
			for (int j = 0; j != n_diagonals_; ++j) {
				if (abs(matrix[j][i] - m.matrix[j][i]) > eps) {
					return false;
				}
			}
		}

		for (int i = 0; i != 2 * n_boundary_rows_; ++i) {
			for (int j = 0; j != n_boundary_elements_; ++j) {
				if (abs(boundary_rows[i][j] - m.boundary_rows[i][j]) > eps) {
					return false;
				}
			}
		}

		return true;
	}
	else {
		return false;
	}
}


TriDiagonal TriDiagonal::operator*(const double scalar) {

	TriDiagonal result(*this);

	scalar_multiply_matrix<TriDiagonal>(scalar, result);

	return result;

}


TriDiagonal TriDiagonal::operator*=(const double scalar) {

	scalar_multiply_matrix<TriDiagonal>(scalar, *this);

	return *this;

}


TriDiagonal TriDiagonal::operator+(const TriDiagonal& rhs) {

	TriDiagonal result(*this);

	matrix_addition<TriDiagonal>(*this, rhs, result);

	return result;

}


TriDiagonal TriDiagonal::operator+=(const TriDiagonal& rhs) {

	matrix_addition<TriDiagonal>(*this, rhs, *this);

	return *this;

}


TriDiagonal TriDiagonal::operator-(const TriDiagonal& rhs) {

	return *this + (-1.0) * rhs;

}


TriDiagonal TriDiagonal::operator-=(const TriDiagonal& rhs) {
	
	*this += (-1.0) * rhs;

	return *this;

}


PentaDiagonal PentaDiagonal::operator*(const double scalar) {

	PentaDiagonal result(*this);

	scalar_multiply_matrix<PentaDiagonal>(scalar, result);

	return result;

}


PentaDiagonal PentaDiagonal::operator*=(const double scalar) {

	scalar_multiply_matrix<PentaDiagonal>(scalar, *this);

	return *this;

}


PentaDiagonal PentaDiagonal::operator+(const PentaDiagonal& rhs) {

	PentaDiagonal result(*this);

	matrix_addition<PentaDiagonal>(*this, rhs, result);

	return result;

}


PentaDiagonal PentaDiagonal::operator+=(const PentaDiagonal& rhs) {

	matrix_addition<PentaDiagonal>(*this, rhs, *this);

	return *this;

}


PentaDiagonal PentaDiagonal::operator-(const PentaDiagonal& rhs) {

	return *this + (-1.0) * rhs;

}


PentaDiagonal PentaDiagonal::operator-=(const PentaDiagonal& rhs) {

	*this += (-1.0) * rhs;

	return *this;

}


// Add scalar to main diagonal.
void BandDiagonal::add_diagonal(const double scalar) {

	for (int i = 0; i != order_; ++i) {
		matrix[bandwidth_][i] += scalar;
	}

	for (int i = 0; i != n_boundary_rows_; ++i) {
		boundary_rows[i][i] += scalar;

		int br_index = (2 * n_boundary_rows_ - 1) - i;
		int be_index = (n_boundary_elements_ - 1) - i;

		boundary_rows[br_index][be_index] += scalar;
	}

}


// Add vector to main diagonal.
void BandDiagonal::add_diagonal(const std::vector<double>& diagonal) {

	for (int i = 0; i != order_; ++i) {
		matrix[bandwidth_][i] += diagonal[i];
	}

	for (int i = 0; i != n_boundary_rows_; ++i) {
		boundary_rows[i][i] += diagonal[i];

		int index1 = (n_boundary_elements_ - 1) - i;
		int index2 = (order_ - 1) - i;

		boundary_rows[index1][index1] += diagonal[index2];
	}

}


// Multiply each element with scalar.
void BandDiagonal::scalar_prod(const double scalar) {

	for (int i = 0; i != n_diagonals_; ++i) {
		for (int j = 0; j != order_; ++j) {
			matrix[i][j] *= scalar;
		}
	}

	for (int i = 0; i != 2 * n_boundary_rows_; ++i) {
		for (int j = 0; j != n_boundary_elements_; ++j) {
			boundary_rows[i][j] *= scalar;
		}
	}

}


// Add two identical matrices.
BandDiagonal BandDiagonal::add_matrix(BandDiagonal mat) {

	for (int i = 0; i != n_diagonals_; ++i) {
		for (int j = 0; j != order_; ++j) {
			mat.matrix[i][j] += matrix[i][j];
		}
	}

	for (int i = 0; i != 2 * n_boundary_rows_; ++i) {
		for (int j = 0; j != n_boundary_elements_; ++j) {
			mat.boundary_rows[i][j] += boundary_rows[i][j];
		}
	}

	return mat;

}


std::vector<double> BandDiagonal::mat_vec_prod(const std::vector<double>& column) {

	std::vector<double> result(order_, 0.0);

	int mr_lower_idx = 0;
	int mr_upper_idx = 0;
	int br_lower_idx = 0;
	int br_upper_idx = 0;

	int be_upper_idx = 0;
	int column_idx = 0;

	// Boundary rows.
	for (int i = 0; i != n_boundary_rows_; ++i) {
		for (int j = i; j != n_boundary_elements_; ++j) {

			mr_lower_idx = i;
			mr_upper_idx = (order_ - 1) - mr_lower_idx;

			br_lower_idx = i;
			br_upper_idx = (2 * n_boundary_rows_ - 1) - br_lower_idx;

			// Lower boundary row.
			result[mr_lower_idx] += boundary_rows[br_lower_idx][j] * column[j];

			be_upper_idx = (n_boundary_elements_ - 1) - j;
			column_idx = (order_ - 1) - j;

			// Upper boundary row.
			result[mr_upper_idx] += boundary_rows[br_upper_idx][be_upper_idx] * column[column_idx];

		}
	}

	// Interior rows.
	for (int i = n_boundary_rows_; i != order_ - n_boundary_rows_; ++i) {
		for (int j = 0; j != n_diagonals_; ++j) {
			result[i] += matrix[j][i] * column[(i - n_boundary_rows_) + j];
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
		boundary_rows_tmp[br_upper_idx][be_upper_idx_tmp] -= upper * matrix[me_upper_idx_tmp][mr_upper_idx];

	}

	// Adjust RHS column vector.
	const int cr_lower_idx = br_lower_idx;
	const int cr_upper_idx = (order_ - 1) - cr_lower_idx;
	column[cr_lower_idx] -= lower * column[mr_lower_idx];
	column[cr_upper_idx] -= upper * column[mr_upper_idx];

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

	int threshold = boundary_row_idx + (bandwidth_ + 1);

	if (threshold > n_boundary_elements_) {
		threshold = n_boundary_elements_;
	}

	// Initialize boundary row.
	for (int i = 0; i != n_diagonals_; ++i) {
		matrix[i][mr_lower_idx] = 0.0;
		matrix[i][mr_upper_idx] = 0.0;
	}

	for (int i = 0; i != threshold; ++i) {

		int me_lower_idx = (bandwidth_ - boundary_row_idx) + i;
		int me_upper_idx = (n_diagonals_ - 1) - me_lower_idx;

		int be_lower_idx = i;
		int be_upper_idx = (n_boundary_elements_ - 1) - be_lower_idx;

		matrix[me_lower_idx][mr_lower_idx] = boundary_rows_tmp[br_lower_idx][be_lower_idx];

		matrix[me_upper_idx][mr_upper_idx] = boundary_rows_tmp[br_upper_idx][be_upper_idx];

	}

}


// Adjust matrix rows at boundary using Gauss elimination.
void TriDiagonal::adjust_boundary(std::vector<double>& column) {

	if (n_boundary_rows_ != 1) {

		throw std::invalid_argument("Number of boundary rows should be 1.");

	}

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
	else {

		throw std::invalid_argument("Number of boundary row elements should be larger than 1 and smaller than 6.");

	}

}


// Adjust matrix rows at boundary using Gauss elimination.
void PentaDiagonal::adjust_boundary(std::vector<double>& column) {

	if (n_boundary_rows_ != 2) {

		throw std::invalid_argument("Number of boundary rows should be 2.");

	}

	if (n_boundary_elements_ >= 2 && n_boundary_elements_ <= 4) {

		boundary_rows_tmp = boundary_rows;
		overwrite_bounary_row(1);
		overwrite_bounary_row(0);

	}
	else if (n_boundary_elements_ == 5) {

		boundary_rows_tmp = boundary_rows;
		gauss_elimination(1, 4, 2, column);
		overwrite_bounary_row(1);
		gauss_elimination(0, 3, 1, column);
		overwrite_bounary_row(0);

	}
	else if (n_boundary_elements_ == 6) {

		boundary_rows_tmp = boundary_rows;
		gauss_elimination(1, 5, 3, column);
		gauss_elimination(1, 4, 2, column);
		overwrite_bounary_row(1);
		gauss_elimination(0, 4, 2, column);
		gauss_elimination(0, 3, 1, column);
		overwrite_bounary_row(0);

	}
	else if (n_boundary_elements_ == 7) {

		boundary_rows_tmp = boundary_rows;
		overwrite_bounary_row(1);
		overwrite_bounary_row(0);

	}
	else {

		throw std::invalid_argument("Number of boundary row elements should be larger than 1 and smaller than 7.");

	}

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


TriDiagonal operator*(const double scalar, TriDiagonal rhs) {

	return rhs * scalar;

}


PentaDiagonal operator*(const double scalar, PentaDiagonal rhs) {

	return rhs * scalar;

}


// TODO: Should the function be void?
template<class T>
void scalar_multiply_matrix(const double scalar, T& matrix) {

	int i_initial = matrix.n_boundary_rows();
	int i_final = matrix.order() - matrix.n_boundary_rows();
	int j_initial = 0;
	int j_final = matrix.n_diagonals();

	for (int i = i_initial; i != i_final; ++i) {
		for (int j = j_initial; j != j_final; ++j) {
			matrix.matrix[j][i] *= scalar;
		}
	}

	i_initial = 0;
	i_final = 2 * matrix.n_boundary_rows();
	j_initial = 0;
	j_final = matrix.n_boundary_elements();

	for (int i = i_initial; i != i_final; ++i) {
		for (int j = j_initial; j != j_final; ++j) {
			matrix.boundary_rows[i][j] *= scalar;
		}
	}

}

// TODO: Should the function be void?
template<class T>
void matrix_addition(const T& matrix1, const T& matrix2, T& result) {

	int i_initial = result.n_boundary_rows();
	int i_final = result.order() - result.n_boundary_rows();
	int j_initial = 0;
	int j_final = result.n_diagonals();

	for (int i = i_initial; i != i_final; ++i) {
		for (int j = j_initial; j != j_final; ++j) {
			result.matrix[j][i] = matrix1.matrix[j][i] + matrix2.matrix[j][i];
		}
	}

	i_initial = 0;
	i_final = 2 * result.n_boundary_rows();
	j_initial = 0;
	j_final = result.n_boundary_elements();

	for (int i = i_initial; i != i_final; ++i) {
		for (int j = j_initial; j != j_final; ++j) {
			result.boundary_rows[i][j] = matrix1.boundary_rows[i][j] + matrix2.boundary_rows[i][j];
		}
	}

}
