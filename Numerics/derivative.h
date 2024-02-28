#pragma once

#include <vector>

#include "band_diagonal_matrix.h"
#include "utility.h"


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
class MixedDerivative {

private:

	T1 d1dx1;
	T2 d1dy1;
	std::vector<double> prefactors;

public:

	MixedDerivative(
		T1& d1dx1_, 
		T2& d1dy1_) {

		d1dx1 = d1dx1_;
		d1dy1 = d1dy1_;

		std::vector<double> tmp(d1dx1.order() * d1dy1.order(), 1.0);
		prefactors = tmp;

	}

	std::vector<double> d2dxdy(
		std::vector<double> func) {

		const int n_x = d1dx1.order();
		const int n_y = d1dy1.order();

		std::vector<double> vec_x(n_x, 0.0);
		std::vector<double> vec_y(n_y, 0.0);

		// Evaluate partial derivative wrt y.
		func = action_2d(n_y, n_x, 2, false, d1dy1, func);

		// Evaluate partial derivative wrt x.
		func = action_2d(n_x, n_y, 1, false, d1dx1, func);

		// Multiply prefactors.
		for (int i = 0; i != n_x * n_y; ++i) {
			func[i] *= prefactors[i];
		}

		return func;

	}

};
