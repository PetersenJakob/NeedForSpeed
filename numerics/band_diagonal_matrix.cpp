#include "band_diagonal_matrix.h"

// Constructor.
BandDiagonal::BandDiagonal(
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

	for (int i = 0; i != matrix.order(); ++i) {
		for (int j = 0; j != matrix.n_diag(); ++j) {
			std::cout << std::setw(14) << matrix.m[j][i];
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

}
