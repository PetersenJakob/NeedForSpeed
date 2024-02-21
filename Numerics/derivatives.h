#pragma once

#include <vector>

struct element {

	// TODO: int, unsigned int, or size_t?
	
	// Row index; first dimension.
	int idx_row;
	// Column index; second dimension.
	int idx_col;
	// Value of element.
	double value;

};

// Finite difference approximation of first order differential operator.
namespace d1dx1 {

	// Central difference; 2nd order accuracy, 3-point stensil.
	std::vector<element> central_2o_3p(int size, double dx);

	// Central difference; 4th order accuracy, 5-point stensil.
	std::vector<element> central_4o_5p(int size, double dx);

}

std::vector<double> mat_vec_product(const int stensil, const std::vector<element>& mat, const std::vector<double>& vec);
