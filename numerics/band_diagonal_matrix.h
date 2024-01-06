#pragma once

#include <iomanip>
#include <iostream>
#include <vector>

// Band-diagonal matrix in compact form. The matrix is assumed to be square.
class BandDiagonal {

	// Order of matrix (number of elements along main diagonal).
	int order_;
	// Number of sub-diagonals (below main diagonal).
	int n_sub_;
	// Number of super-diagonals (above main diagonal).
	int n_super_;
	// Total number of diagonals (including main diagonal).
	int n_diag_;

public:

	// Matrix in compact form.
	std::vector<std::vector<double>> m;

	// Constructor.
	BandDiagonal(
		const int _order,
		const int _n_sub,
		const int _n_super);

	const int order() {
		return order_;
	}

	const int n_sub() {
		return n_sub_;
	}

	const int n_super() {
		return n_super_;
	}

	const int n_diag() {
		return n_diag_;
	}

};

class TriDiagonal : public BandDiagonal {

public:

	TriDiagonal(const int _order) : BandDiagonal(_order, 1, 1) {}

};

class PentaDiagonal : public BandDiagonal {

public:

	PentaDiagonal(const int _order) : BandDiagonal(_order, 2, 2) {}

};

void print_matrix(BandDiagonal matrix);
