#pragma once

#include <vector>


// Band-diagonal matrix stored in compact form.
// Note: The lower and upper bandwidths are assumed to be identical.
class BandDiagonal {

protected:

	// Matrix order: Number of elements along main diagonal.
	int order_;
	// Bandwidth: Number of sub-diagonals or super-diagonals.
	int bandwidth_;
	// Number of diagonals.
	int n_diagonals_;
	// Number of boundary rows (at each boundary).
	int n_boundary_rows_;
	// Number of non-zero elements along each boundary row.
	int n_boundary_elements_;

public:

	// Band-diagonal matrix in compact form.
	std::vector<std::vector<double>> matrix;

	// Boundary rows of band-diagonal matrix.
	std::vector<std::vector<double>> boundary_rows;

	std::vector<std::vector<double>> boundary_rows_tmp;

	BandDiagonal(
		const int _order,
		const int _bandwidth,
		const int _n_boundary_rows,
		const int _n_boundary_elements);


	bool operator == (const BandDiagonal& m)
	{
		const double eps = 1.0e-12;

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


	int order() const {
		return order_;
	}

	int bandwidth() const {
		return bandwidth_;
	}

	int n_diagonals() const {
		return n_diagonals_;
	}

	int n_boundary_rows() const {
		return n_boundary_rows_;
	}

	int n_boundary_elements() const {
		return n_boundary_elements_;
	}

	// Add scalar to main diagonal.
	void add_diagonal(const double scalar);

	// Add vector to main diagonal.
	void add_diagonal(const std::vector<double>& diagonal);

	// Multiply each element with scalar.
	void scalar_prod(const double scalar);

	// Add two identical matrices.
	BandDiagonal add_matrix(BandDiagonal mat);

	// Matrix-vector product.
	std::vector<double> mat_vec_prod(const std::vector<double>& column);

	// ...
	void gauss_elimination(
		const int boundary_row_idx,
		const int boundary_element_idx,
		const int matrix_row_idx,
		std::vector<double>& column);

	// ...
	void overwrite_bounary_row(const int boundary_row_idx);

	// ...
	virtual void adjust_boundary(std::vector<double>& column) {}
	
};

// Tri-diagonal matrix stored in compact form.
class TriDiagonal : public BandDiagonal {

public:

	TriDiagonal(
		const int _order,
		const int _n_boundary_rows = 1,
		const int _n_boundary_elements = 3) : BandDiagonal(_order, 1, _n_boundary_rows, _n_boundary_elements) {}

	void adjust_boundary(std::vector<double>& column);

};

// Penta-diagonal matrix stored in compact form.
class PentaDiagonal : public BandDiagonal {

public:

	PentaDiagonal(
		const int _order,
		const int _n_boundary_rows = 2,
		const int _n_boundary_elements = 3) : BandDiagonal(_order, 2, _n_boundary_rows, _n_boundary_elements) {}

	void adjust_boundary(std::vector<double>& column);

};

void print_matrix(BandDiagonal matrix);
