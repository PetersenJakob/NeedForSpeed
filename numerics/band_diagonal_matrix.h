#pragma once

#include <vector>


// Band-diagonal matrix stored in compact form.
class BandDiagonal {

	// Matrix order: Number of elements along main diagonal.
	int order_;
	// Lower bandwidth: Number of non-zero sub-diagonals.
	int lower_bandwidth_;
	// Upper bandwidth: Number of non-zero super-diagonals.
	int upper_bandwidth_;
	// Bandwidth: Number of non-zero diagonals, including main diagonal.
	int bandwidth_;

	// Number of boundary rows.
	int n_boundary_rows_;
	// Number of non-zero elements along boundary rows.
	int n_boundary_elements_;

public:

	BandDiagonal(
		const int _order,
		const int _lower_bandwidth,
		const int _upper_bandwidth,
		const int _n_boundary_rows,
		const int _n_boundary_elements);

	const int order() {
		return order_;
	}

	const int lower_bandwidth() {
		return lower_bandwidth_;
	}

	const int upper_bandwidth() {
		return upper_bandwidth_;
	}

	const int bandwidth() {
		return bandwidth_;
	}

	const int n_boundary_rows() {
		return n_boundary_rows_;
	}

	const int n_boundary_elements() {
		return n_boundary_elements_;
	}

	// Band-diagonal matrix in compact form.
	std::vector<std::vector<double>> matrix;

	// Boundary rows of band-diagonal matrix (ordered according to row indices).
	std::vector<std::vector<double>> boundary_rows;

};

// Tri-diagonal matrix stored in compact form.
class TriDiagonal : public BandDiagonal {

public:

	TriDiagonal(
		const int _order,
		const int _n_boundary_rows = 1,
		const int _n_boundary_elements = 3) : BandDiagonal(_order, 1, 1, _n_boundary_rows, _n_boundary_elements) {}

};

// Penta-diagonal matrix stored in compact form.
class PentaDiagonal : public BandDiagonal {

public:

	PentaDiagonal(
		const int _order,
		const int _n_boundary_rows = 2,
		const int _n_boundary_elements = 5) : BandDiagonal(_order, 2, 2, _n_boundary_rows, _n_boundary_elements) {}

};

void print_matrix(BandDiagonal matrix);
