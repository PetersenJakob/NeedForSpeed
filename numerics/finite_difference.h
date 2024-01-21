#pragma once

#include "band_diagonal_matrix.h"


// Finite difference representations of first order derivative operator.
namespace d1dx1 {

	// Central difference; 2nd order accuracy. Boundary; 1st order accuracy.
	TriDiagonal c2b1(const int order, const double dx);

	// Central difference; 2nd order accuracy. Boundary; 2nd order accuracy.
	TriDiagonal c2b2(const int order, const double dx);

	// Central difference; 4th order accuracy. Boundary; 4th order accuracy.
	PentaDiagonal c4b4(const int order, const double dx);

}

// Finite difference representations of second order derivative operator.
namespace d2dx2 {

	// Central difference; 2nd order accuracy. Boundary; 1st order accuracy.
	TriDiagonal c2b1(const int order, const double dx);

	// Central difference; 2nd order accuracy. Boundary; 2nd order accuracy.
	TriDiagonal c2b2(const int order, const double dx);

	// Central difference; 4th order accuracy. Boundary; 4th order accuracy.
	PentaDiagonal c4b4(const int order, const double dx);

}

// Adjusting finite difference representation at boundary.
template <class T>
void boundary(const int row_index, const double dx, const std::vector<double>& coef, T& matrix);
