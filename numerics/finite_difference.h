#pragma once

#include "band_diagonal_matrix.h"

namespace d1dx1 {

	// Central difference; 2nd order accuracy, 3-point stensil.
	TriDiagonal central_tri(const int order, const double dx);

	// Central differnce; 4th order accuracy, 5-point stensil.
	PentaDiagonal central_penta(const int order, const double dx);

	// Forward difference; 1st order accuracy, 2-point stensil.
	TriDiagonal forward_1o(const int order, const double dx);

	// Forward difference; 2nd order accuracy, 3-point stensil.
	PentaDiagonal forward_2o(const int order, const double dx);

	// Backward difference; 1st order accuracy, 2-point stensil.
	TriDiagonal backward_1o(const int order, const double dx);

	// Backward difference; 2nd order accuracy, 3-point stensil.
	PentaDiagonal backward_2o(const int order, const double dx);

}