#pragma once

#include <vector>

#include "band_diagonal_matrix.h"
#include "coefficients.h"
#include "utility.h"


// Finite difference representation of first order derivative operator.
namespace d1dx1Template {

	// Finite difference representation on uniform grid.
	namespace uniform {

		// Interior: Central difference, 2nd order accuracy.
		// Boundary 1st row: Forward difference, 1st order accurary.
		template<typename Tnumber>
		TriDiagonalTemplate<Tnumber> c2b1(const std::vector<Tnumber>& grid);

		// Interior: Central difference, 2nd order accuracy.
		// Boundary 1st row: Forward difference, 2nd order accurary.
		template<typename Tnumber>
		TriDiagonalTemplate<Tnumber> c2b2(const std::vector<Tnumber>& grid);

		// Interior: Central difference, 4th order accuracy.
		// Boundary 1st row: Forward difference, 2nd order accurary.
		// Boundary 2nd row: Central difference, 2nd order accurary.
		template<typename Tnumber>
		PentaDiagonalTemplate<Tnumber> c4b2(const std::vector<Tnumber>& grid);

		// Interior: Central difference, 4th order accuracy.
		// Boundary 1st row: Forward difference, 4th order accurary.
		// Boundary 2nd row: Forward difference, 4th order accurary.
		template<typename Tnumber>
		PentaDiagonalTemplate<Tnumber> c4b4(const std::vector<Tnumber>& grid);

	}
}


// First order derivative operator. 
// Central difference; 2nd order accuracy. Boundary; 1st order accuracy.
template<typename Tnumber>
TriDiagonalTemplate<Tnumber> d1dx1Template::uniform::c2b1(const std::vector<Tnumber>& grid) {

	using MatrixType = TriDiagonalTemplate<Tnumber>;

	// TODO: If you choose two boundary rows with f1+f1 and b1+b1, what will the L2 function norm be?

	const std::size_t order = grid.size();
	// TODO: Should dx be of type Tnumber?
	const Tnumber dx = grid[1] - grid[0];

	MatrixType matrix = 
		setupTemplate<MatrixType, Tnumber>(order, coef_x1Template::uniform::c2<Tnumber>(dx), 1, 2);
	boundaryTemplate<MatrixType, Tnumber>(0, "lower", coef_x1Template::uniform::f1<Tnumber>(dx), matrix);
	boundaryTemplate<MatrixType, Tnumber>(0, "upper", coef_x1Template::uniform::b1<Tnumber>(dx), matrix);

	return matrix;

}


// First order derivative operator. 
// Central difference; 2nd order accuracy. Boundary; 2nd order accuracy.
template<typename Tnumber>
TriDiagonalTemplate<Tnumber> d1dx1Template::uniform::c2b2(const std::vector<Tnumber>& grid) {

	using MatrixType = TriDiagonalTemplate<Tnumber>;

	const std::size_t order = grid.size();
	const Tnumber dx = grid[1] - grid[0];

	MatrixType matrix = 
		setup<MatrixType, Tnumber>(order, coef_x1Template::uniform::c2<Tnumber>(dx), 1, 3);
	boundaryTemplate<MatrixType, Tnumber>(0, "lower", coef_x1Template::uniform::f2<Tnumber>(dx), matrix);
	boundaryTemplate<MatrixType, Tnumber>(0, "upper", coef_x1Template::uniform::b2<Tnumber>(dx), matrix);

	return matrix;

}


// First order derivative operator. 
// Central difference; 4th order accuracy. Boundary; 2nd order accuracy.
template<typename Tnumber>
PentaDiagonalTemplate<Tnumber> d1dx1Template::uniform::c4b2(const std::vector<Tnumber>& grid) {

	using MatrixType = PentaDiagonalTemplate<Tnumber>;

	const std::size_t order = grid.size();
	const Tnumber dx = grid[1] - grid[0];

	MatrixType matrix =
		setupTemplate<MatrixType, Tnumber>(order, coef_x1Template::uniform::c4<Tnumber>(dx), 2, 3);
	boundaryTemplate<MatrixType, Tnumber>(0, "lower", coef_x1Template::uniform::f2<Tnumber>(dx), matrix);
	boundaryTemplate<MatrixType, Tnumber>(1, "lower", coef_x1Template::uniform::f2<Tnumber>(dx), matrix);
	boundaryTemplate<MatrixType, Tnumber>(1, "upper", coef_x1Template::uniform::b2<Tnumber>(dx), matrix);
	boundaryTemplate<MatrixType, Tnumber>(0, "upper", coef_x1Template::uniform::b2<Tnumber>(dx), matrix);

	return matrix;

}


// First order derivative operator. 
// Central difference; 4th order accuracy. Boundary; 4th order accuracy.
template<typename Tnumber>
PentaDiagonalTemplate<Tnumber> d1dx1Template::uniform::c4b4(const std::vector<Tnumber>& grid) {

	using MatrixType = PentaDiagonalTemplate<Tnumber>;

	const std::size_t order = grid.size();
	const Tnumber dx = grid[1] - grid[0];

	MatrixType matrix = setup<MatrixType, Tnumber>(order, coef_x1Template::uniform::c4<Tnumber>(dx), 2, 5);
	boundaryTemplate<MatrixType, Tnumber>(0, "lower", coef_x1Template::uniform::f4<Tnumber>(dx), matrix);
	boundaryTemplate<MatrixType, Tnumber>(1, "lower", coef_x1Template::uniform::f4<Tnumber>(dx), matrix);
	boundaryTemplate<MatrixType, Tnumber>(1, "upper", coef_x1Template::uniform::b4<Tnumber>(dx), matrix);
	boundaryTemplate<MatrixType, Tnumber>(0, "upper", coef_x1Template::uniform::b4<Tnumber>(dx), matrix);

	return matrix;

}


// ###############################################################################


// Finite difference representation of first order derivative operator.
namespace d1dx1 {

	// Finite difference representation on uniform grid.
	namespace uniform {

		// Interior: Central difference, 2nd order accuracy.
		// Boundary 1st row: Forward difference, 1st order accurary.
		TriDiagonal c2b1(const std::vector<double>& grid);

		// Interior: Central difference, 2nd order accuracy.
		// Boundary 1st row: Forward difference, 2nd order accurary.
		TriDiagonal c2b2(const std::vector<double>& grid);

		// Interior: Central difference, 4th order accuracy.
		// Boundary 1st row: Forward difference, 2nd order accurary.
		// Boundary 2nd row: Central difference, 2nd order accurary.
		PentaDiagonal c4b2(const std::vector<double>& grid);

		// Interior: Central difference, 4th order accuracy.
		// Boundary 1st row: Forward difference, 4th order accurary.
		// Boundary 2nd row: Forward difference, 4th order accurary.
		PentaDiagonal c4b4(const std::vector<double>& grid);

	}

	// Finite difference representation on non-uniform grid.
	namespace nonuniform {

		// Interior: Central difference, 2nd order accuracy.
		// Boundary 1st row: Forward difference, 1st order accurary.
		TriDiagonal c2b1(const std::vector<double>& grid);

		// Interior: Central difference, 2nd order accuracy.
		// Boundary 1st row: Forward difference, 2nd order accurary.
		TriDiagonal c2b2(const std::vector<double>& grid);

		// Interior: Central difference, 4th order accuracy.
		// Boundary 1st row: Forward difference, 2nd order accurary.
		// Boundary 2nd row: Central difference, 2nd order accurary.
		PentaDiagonal c4b2(const std::vector<double>& grid);

	}

}


// Finite difference representation of second order derivative operator.
namespace d2dx2 {

	// Finite difference representation on uniform grid.
	namespace uniform {

		// Interior: Central difference, 2nd order accuracy.
		// Boundary 1st row: Neumann boundary condition, d2dx2 = 0. TODO: Is this Neumann, or only the case for B with fist order derivatives.
		TriDiagonal c2b0(const std::vector<double>& grid);

		// Interior: Central difference, 2nd order accuracy.
		// Boundary 1st row: Forward difference, 1st order accurary.
		TriDiagonal c2b1(const std::vector<double>& grid);

		// Interior: Central difference, 2nd order accuracy.
		// Boundary 1st row: Forward difference, 2nd order accurary.
		TriDiagonal c2b2(const std::vector<double>& grid);

		// Interior: Central difference, 4th order accuracy.
		// Boundary 1st row: Neumann boundary condition, d2dx2 = 0.
		// Boundary 2nd row: Central difference, 2nd order accurary.
		PentaDiagonal c4b0(const std::vector<double>& grid);

		// Interior: Central difference, 4th order accuracy.
		// Boundary 1st row: Forward difference, 2nd order accurary.
		// Boundary 2nd row: Central difference, 2nd order accurary.
		PentaDiagonal c4b2(const std::vector<double>& grid);

		// Interior: Central difference, 4th order accuracy.
		// Boundary 1st row: Forward difference, 4th order accurary.
		// Boundary 2nd row: Forward difference, 4th order accurary.
		PentaDiagonal c4b4(const std::vector<double>& grid);

	}

	// Finite difference representation on non-uniform grid.
	namespace nonuniform {

		// Interior: Central difference, 2nd order accuracy.
		// Boundary 1st row: Neumann boundary condition, d2dx2 = 0.
		TriDiagonal c2b0(const std::vector<double>& grid);

		// Interior: Central difference, 2nd order accuracy.
		// Boundary 1st row: Forward difference, 1st order accurary.
		TriDiagonal c2b1(const std::vector<double>& grid);

		// Interior: Central difference, 4th order accuracy.
		// Boundary 1st row: Neumann boundary condition, d2dx2 = 0.
		// Boundary 2nd row: Central difference, 2nd order accurary.
		PentaDiagonal c4b0(const std::vector<double>& grid);

	}

}


// TODO: Class method for evaluating second order mixed derivative operator.
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

	void set_prefactors(const double scalar) {
		for (int i = 0; i != prefactors.size(); ++i) {
			prefactors[i] = scalar;
		}
	}

	void set_prefactors(
		const std::vector<double>& coef_x,
		const std::vector<double>& coef_y) {
		int index = 0;
		for (int i = 0; i != coef_x.size(); ++i) {
			for (int j = 0; j != coef_y.size(); ++j) {
				prefactors[i] = coef_x[i] * coef_y[j];
				++index;
			}
		}
	}

	void set_prefactors(
		const std::vector<double>& factors) {
		int index = 0;
		for (int i = 0; i != factors.size(); ++i) {
			prefactors[i] = factors[i];
				++index;
		}
	}

	std::vector<double> d2dxdy(
		std::vector<double> func) {

		const int n_x = d1dx1.order();
		const int n_y = d1dy1.order();

		std::vector<double> vec_x(n_x, 0.0);
		std::vector<double> vec_y(n_y, 0.0);

		// Evaluate first order partial derivative wrt y.
		func = action_2d(n_y, n_x, 2, false, d1dy1, func);

		// Evaluate first order partial derivative wrt x.
		func = action_2d(n_x, n_y, 1, false, d1dx1, func);

		// Multiply prefactors.
		for (int i = 0; i != n_x * n_y; ++i) {
			func[i] *= prefactors[i];
		}

		return func;

	}

};
