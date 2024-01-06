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
		const int _n_super) {

		order_ = _order;
		n_sub_ = _n_sub;
		n_super_ = _n_super;
		n_diag_ = 1 + n_sub_ + n_super_;

		// First index: Diagonal (sub -> main -> super).
		// Second index: Row index.
		std::vector<double> diagonal(order_, 0.0);
		std::vector<std::vector<double>> matrix(n_diag_, diagonal);

#if false
		// First index: Row index.
		// Second index: Diagonal (sub -> main -> super).
		std::vector<double> row(n_diag_, 0.0);
		std::vector<std::vector<double>> matrix(order_, row);
#endif

		m = matrix;

	}

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

// TODO: Why not const parameter allowed?
void print_matrix(BandDiagonal& matrix) {

	std::cout << std::scientific << std::setprecision(10);

	std::cout 
		<< "Band-diagonal matrix (n_sub = " 
		<< matrix.n_sub() 
		<< ", 1, n_super = " 
		<< matrix.n_super() 
		<< ") of order " 
		<< matrix.order() << ":" << std::endl;

	for (int i = 0; i != matrix.order(); ++i) {
		for (int j = 0; j != matrix.n_diag(); ++j) {
			std::cout << std::setw(20) << matrix.m[j][i];
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

}
