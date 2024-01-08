#include "band_diagonal_matrix.h"

// Constructor.
BandDiagonal::BandDiagonal(
	const int _order,
	const int _n_sub,
	const int _n_super,
	const int _n_bound_rows,
	const int _n_bound_elements) {

	order_ = _order;
	n_sub_ = _n_sub;
	n_super_ = _n_super;
	n_diag_ = 1 + n_sub_ + n_super_;

	n_bound_rows_ = _n_bound_rows;
	n_bound_elements_ = _n_bound_elements;

	// First index: Diagonal (sub -> main -> super).
	// Second index: Row index.
	std::vector<double> diagonal(order_, 0.0);
	std::vector<std::vector<double>> vec1(n_diag_, diagonal);
	m = vec1;

	std::vector<double> row(n_bound_elements_);
	std::vector<std::vector<double>> vec2(2 * n_bound_rows_, row);
	b = vec2;

}

// TODO: Why not const or reference allowed?
void print_matrix(BandDiagonal matrix) {

	std::cout << std::scientific << std::setprecision(5);

	std::cout
		<< "Band-diagonal matrix (n_sub = "
		<< matrix.n_sub()
		<< ", 1, n_super = "
		<< matrix.n_super()
		<< ") of order "
		<< matrix.order() << ":" << std::endl;

	std::cout << "Boundary rows" << std::endl;
	for (int i = 0; i != matrix.n_bound_rows(); ++i) {
		for (int j = 0; j != matrix.n_bound_elements(); ++j) {
			std::cout << std::setw(14) << matrix.b[i][j];
		}
		std::cout << std::endl;
	}

	std::cout << "Matrix" << std::endl;
	for (int i = 0; i != matrix.order(); ++i) {
		for (int j = 0; j != matrix.n_diag(); ++j) {
			std::cout << std::setw(14) << matrix.m[j][i];
		}
		std::cout << std::endl;
	}

	std::cout << "Boundary rows" << std::endl;
	for (int i = matrix.n_bound_rows(); i != 2 * matrix.n_bound_rows(); ++i) {
		for (int j = 0; j != matrix.n_bound_elements(); ++j) {
			std::cout << std::setw(14) << matrix.b[i][j];
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

}
