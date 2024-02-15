#pragma once

#include "band_diagonal_matrix.h"


// Matrix representation of identity operator.
namespace identity {

	TriDiagonal tri(const int order, const int n_boundary_elements = 2);

	PentaDiagonal penta(const int order, const int n_boundary_elements = 3);

}


// Finite difference representation of first order derivative operator.
namespace d1dx1 {

	// Finite difference representation on uniform grid.
	namespace uniform {

		// Central difference; 2nd order accuracy. Boundary; 1st order accuracy.
		TriDiagonal c2b1(const int order, const double dx);

		// Central difference; 2nd order accuracy. Boundary; 2nd order accuracy.
		TriDiagonal c2b2(const int order, const double dx);


		// TODO: Test convergence rate
		// Central difference; 2nd order accuracy. Boundary; 3rd order accuracy.
		TriDiagonal c2b3(const int order, const double dx);


		// Central difference; 4th order accuracy. Boundary; 2nd order accuracy.
		PentaDiagonal c4b2(const int order, const double dx);

		// Central difference; 4th order accuracy. Boundary; 4th order accuracy.
		PentaDiagonal c4b4(const int order, const double dx);

	}

	// Finite difference representation on non-uniform grid.
	namespace nonuniform {

		// Central difference; 2nd order accuracy. Boundary; 1st order accuracy.
		TriDiagonal c2b1(const int order, const std::vector<double> grid);

		// Central difference; 2nd order accuracy. Boundary; 2nd order accuracy.
		TriDiagonal c2b2(const int order, const std::vector<double> grid);

		// Central difference; 4th order accuracy. Boundary; 2nd order accuracy.
		PentaDiagonal c4b2(const int order, const std::vector<double> grid);

	}

}


// Finite difference representation of second order derivative operator.
namespace d2dx2 {

	// Finite difference representation on uniform grid.
	namespace uniform {

		// Central difference; 2nd order accuracy. Boundary; d2dx2 = 0.
		TriDiagonal c2b0(const int order, const double dx);

		// Central difference; 2nd order accuracy. Boundary; 1st order accuracy.
		TriDiagonal c2b1(const int order, const double dx);

		// Central difference; 2nd order accuracy. Boundary; 2nd order accuracy.
		TriDiagonal c2b2(const int order, const double dx);


		// Central difference; 4th order accuracy. Boundary; 2nd row c2, 1st row d2dx2 = 0.
		PentaDiagonal c4b0(const int order, const double dx);


		// Central difference; 4th order accuracy. Boundary; 4th order accuracy.
		PentaDiagonal c4b4(const int order, const double dx);

	}

	// Finite difference representation on non-uniform grid.
	namespace nonuniform {

		// Central difference; 2nd order accuracy. Boundary; d2dx2 = 0.
		TriDiagonal c2b0(const int order, const std::vector<double> grid);

		// Central difference; 2nd order accuracy. Boundary; 1st order accuracy.
		TriDiagonal c2b1(const int order, const std::vector<double> grid);

		// Central difference; 4th order accuracy. Boundary; 2nd row c2, 1st row d2dx2 = 0.
		PentaDiagonal c4b0(const int order, const std::vector<double> grid);

	}

}


// Adjusting finite difference representation at boundary.
template <class T>
void boundary(const int row_index, const std::vector<double>& coef, T& matrix);
