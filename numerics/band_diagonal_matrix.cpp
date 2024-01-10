#include <iomanip>
#include <iostream>

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
