#pragma once

#include <vector>

#include "band_diagonal_matrix.h"


// Finite difference representation of first order derivative operator.
namespace d1dx1 {

	// TODO: Which ones to keep?

	// Finite difference representation on uniform grid.
	namespace uniform {

		// Interior: Central difference, 2nd order accuracy. 
		// Boundary: 1st order accuracy.
		TriDiagonal c2b1(const int order, const double dx);

		// Interior: Central difference, 2nd order accuracy. 
		// Boundary: 2nd order accuracy.
		TriDiagonal c2b2(const int order, const double dx);

		// Interior: Central difference, 4th order accuracy.
		// Boundary: 2nd order accuracy.
		PentaDiagonal c4b2(const int order, const double dx);

		// Interior: Central difference, 4th order accuracy. 
		// Boundary: 4th order accuracy.
		PentaDiagonal c4b4(const int order, const double dx);

	}

	// Finite difference representation on non-uniform grid.
	namespace nonuniform {

		// Interior: Central difference, 2nd order accuracy. 
		// Boundary: 1st order accuracy.
		TriDiagonal c2b1(const int order, const std::vector<double> grid);

		// Interior: Central difference, 2nd order accuracy. 
		// Boundary: 2nd order accuracy.
		TriDiagonal c2b2(const int order, const std::vector<double> grid);

		// Interior: Central difference, 4th order accuracy. 
		// Boundary. 2nd order accuracy.
		PentaDiagonal c4b2(const int order, const std::vector<double> grid);

	}

}


// Finite difference representation of second order derivative operator.
namespace d2dx2 {

	// TODO: Which ones to keep?

	// Finite difference representation on uniform grid.
	namespace uniform {

		// Interior: Central difference, 2nd order accuracy. 
		// Boundary; d2dx2 = 0.
		TriDiagonal c2b0(const int order, const double dx);

		// Interior: Central difference, 2nd order accuracy. 
		// Boundary; 1st order accuracy.
		TriDiagonal c2b1(const int order, const double dx);

		// Interior: Central difference, 2nd order accuracy. 
		// Boundary; 2nd order accuracy.
		TriDiagonal c2b2(const int order, const double dx);

		// Interior: Central difference, 4th order accuracy. 
		// Boundary; 2nd row c2, 1st row d2dx2 = 0.
		PentaDiagonal c4b0(const int order, const double dx);

		// Interior: Central difference, 4th order accuracy. 
		// Boundary; 4th order accuracy.
		PentaDiagonal c4b4(const int order, const double dx);

	}

	// Finite difference representation on non-uniform grid.
	namespace nonuniform {

		// Interior: Central difference, 2nd order accuracy. 
		// Boundary; d2dx2 = 0.
		TriDiagonal c2b0(const int order, const std::vector<double> grid);

		// Interior: Central difference, 2nd order accuracy. 
		// Boundary; 1st order accuracy.
		TriDiagonal c2b1(const int order, const std::vector<double> grid);

		// Interior: Central difference, 4th order accuracy. 
		// Boundary; 2nd row c2, 1st row d2dx2 = 0.
		PentaDiagonal c4b0(const int order, const std::vector<double> grid);

	}

}


// Finite difference representation of second order mixed derivative operator.
template <class T1, class T2>
std::vector<std::vector<double>> d2dxdy(
	T1& d1dx1,
	T2& d1dy1,
	std::vector<std::vector<double>> func) {

	const int n_points_x = (int)func.size();
	const int n_points_y = (int)func[0].size();

	std::vector<double> vec_x(n_points_x, 0.0);
	std::vector<double> vec_y(n_points_y, 0.0);

	// Evaluate partial derivative wrt y.
	for (int i = 0; i != n_points_x; ++i) {

		for (int j = 0; j != n_points_y; ++j) {
			vec_y[j] = func[i][j];
		}

		vec_y = d1dy1 * vec_y;

		for (int j = 0; j != n_points_y; ++j) {
			func[i][j] = vec_y[j];
		}

	}

	// Evaluate partial derivative wrt x.
	for (int i = 0; i != n_points_y; ++i) {

		for (int j = 0; j != n_points_x; ++j) {
			vec_x[j] = func[j][i];
		}

		vec_x = d1dx1 * vec_x;

		for (int j = 0; j != n_points_x; ++j) {
			func[j][i] = vec_x[j];
		}

	}

	return func;

}
