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

	BandDiagonal(
		const int _order,
		const int _bandwidth,
		const int _n_boundary_rows,
		const int _n_boundary_elements);


	const int order() {
		return order_;
	}

	const int bandwidth() {
		return bandwidth_;
	}

	const int n_diagonals() {
		return n_diagonals_;
	}

	const int n_boundary_rows() {
		return n_boundary_rows_;
	}

	const int n_boundary_elements() {
		return n_boundary_elements_;
	}

	// Matrix-vector product.
	std::vector<double> mat_vec_prod(const std::vector<double>& column);

	// TODO: Mover from TriDiagonal!
//	void adjust_boundary(std::vector<double>& column);

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

};

void print_matrix(BandDiagonal matrix);
