#pragma once

#include "band_diagonal_matrix.h"


// Finite difference representations of first order derivative operator.
namespace d1dx1 {

	// Central difference; 2nd order accuracy.
	TriDiagonal central_tri(const int order, const double dx);

	// Central difference; 4th order accuracy.
	PentaDiagonal central_penta(const int order, const double dx);

	// Forward difference; 1st order accuracy.
	TriDiagonal f1(const int order, const double dx);

	// Forward difference; 2nd order accuracy.
	PentaDiagonal f2(const int order, const double dx);

	// Backward difference; 1st order accuracy.
	TriDiagonal b1(const int order, const double dx);

	// Backward difference; 2nd order accuracy.
	PentaDiagonal b2(const int order, const double dx);

}

// Finite difference representations of second order derivative operator.
namespace d2dx2 {

	// Central difference; 2nd order accuracy.
	TriDiagonal central_tri(const int order, const double dx);

	// Central difference; 4th order accuracy.
	PentaDiagonal central_penta(const int order, const double dx);

	// Forward difference; 1st order accuracy.
	TriDiagonal f1(const int order, const double dx);

	// Forward difference; 2nd order accuracy.
	PentaDiagonal f2(const int order, const double dx);

	// Backward difference; 1st order accuracy.
	TriDiagonal b1(const int order, const double dx);

	// Backward difference; 2nd order accuracy.
	PentaDiagonal b2(const int order, const double dx);

}
