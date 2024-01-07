//#include <algorithm>
#include <algorithm>
#include <vector>

#include "band_diagonal_matrix.h"
#include "finite_difference.h"


// Finite difference coefficients.
// TODO: Coefficients are always given for a symmetric set of points with respect to the main diagonal.


// Reverse order of coefficients and multiply by scalar.
std::vector<double> reverse_order(std::vector<double> coef, const int scalar) {

	std::reverse(coef.begin(), coef.end());

	for (auto& elem : coef) {
		elem *= scalar;
	}

	return coef;

}

// Coefficients for finite difference representations of first order derivative operator.
namespace coef1 {

    // Central difference; 2nd order accuracy.
    std::vector<double> c2 { 
		 1.0 / 2.0, 
		 0.0, 
		-1.0 / 2.0 
	};

    // Central difference; 4th order accuracy.
    std::vector<double> c4 { 
		-1.0 / 12.0, 
		 2.0 / 3.0, 
		 0.0, 
		-2.0 / 3.0, 
		 1.0 / 12.0 
	};

    // Forward difference; 1st order accuracy.
    std::vector<double> f1 { 
		 1.0, 
		-1.0, 
		 0.0 
	};

    // Forward difference; 2nd order accuracy.
    std::vector<double> f2 { 
		-1.0 / 2.0, 
		 2.0, 
		-3.0 / 2.0, 
		 0.0, 
		 0.0 
	};

	// Forward difference; 3rd order accuracy.
	std::vector<double> f3 {
		 1.0 / 3.0,
		-3.0 / 2.0,
		 3.0,
		-11.0 / 16.0,
		 0.0,
		 0.0,
		 0.0 
	};

	// Forward difference; 4th order accuracy.
	std::vector<double> f4 {
		-1.0 / 4.0,
		 4.0 / 3.0,
		-3.0,
		 4.0,
		-25.0 / 12.0,
		 0.0,
		 0.0,
		 0.0,
		 0.0 
	};

	// Backward difference; 1st order accuracy.
	std::vector<double> b1 = reverse_order(f1, -1.0);

	// Backward difference; 2nd order accuracy.
	std::vector<double> b2 = reverse_order(f2, -1.0);

	// Backward difference; 3rd order accuracy.
	std::vector<double> b3 = reverse_order(f3, -1.0);

	// Backward difference; 4th order accuracy.
	std::vector<double> b4 = reverse_order(f4, -1.0);

}

// Coefficients for finite difference representations of second order derivative operator.
namespace coef2 {

	// Central difference; 2nd order accuracy.
	std::vector<double> c2 { 
		 1.0, 
		-2.0, 
		 1.0 
	};

	// Central difference; 4th order accuracy.
	std::vector<double> c4 { 
		-1.0 / 12.0, 
		 4.0 / 3.0, 
		-5.0 / 2.0, 
		 4.0 / 3.0, 
		-1.0 / 12.0 
	};

	// Forward difference; 1st order accuracy.
	std::vector<double> f1 { 
		 1.0, 
		-2.0, 
		 1.0,
		 0.0,
		 0.0 
	};

	// Forward difference; 2nd order accuracy.
	std::vector<double> f2 { 
		-1.0, 
		 4.0, 
		-5.0, 
		 2.0, 
		 0.0,
		 0.0,
		 0.0 
	};

	// Forward difference; 3rd order accuracy.
	std::vector<double> f3 { 
		 11.0 / 12.0, 
		-14.0 / 3.0, 
		 19.0 / 2.0, 
		-26.0 / 3.0,
		 35.0 / 12.0, 
		 0.0,
		 0.0,
		 0.0,
		 0.0 
	};

	// Forward difference; 4th order accuracy.
	std::vector<double> f4 {
		-5.0 / 6.0,
		 61.0 / 12.0,
		-13.0,
		 107.0 / 6.0,
		-77.0 / 6.0,
		 15.0 / 4.0,
		 0.0,
		 0.0,
		 0.0,
		 0.0,
		 0.0 
	};

	// Backward difference; 1st order accuracy.
	std::vector<double> b1 = reverse_order(f1, 1.0);

	// Backward difference; 2nd order accuracy.
	std::vector<double> b2 = reverse_order(f2, 1.0);

	// Backward difference; 3rd order accuracy.
	std::vector<double> b3 = reverse_order(f3, 1.0);

	// Backward difference; 4th order accuracy.
	std::vector<double> b4 = reverse_order(f4, 1.0);

}

// Setting up finite difference representation of derivative operator on equidistant grid.
template <class T> 
T setup(const int order, const double dx, const std::vector<double>& coef) {

	T matrix(order);

	for (int i = 0; i != coef.size(); ++i) {
		for (int j = 0; j != order; ++j) {
			matrix.m[i][j] = coef[i] / dx;
		}
	}

	return matrix;

}

// Adjusting finite difference representations at boundary.
template <class T>
void boundary(const int row_index, const double dx, const std::vector<double>& coef, T& matrix) {
	
	for (int i = 0; i != coef.size(); ++i) {
		matrix.m[i][row_index] = coef[i] / dx;
	}

}

// Central difference; 2nd order accuracy.
TriDiagonal d1dx1::central_tri(const int order, const double dx) {

	TriDiagonal matrix = setup<TriDiagonal>(order, dx, coef1::c2);

	// Adjust finite difference approximations at boundary.
	boundary(0, dx, coef1::b1, matrix);
	boundary(order - 1, dx, coef1::f1, matrix);

	return matrix;

}

// Central difference; 4th order accuracy.
PentaDiagonal d1dx1::central_penta(const int order, const double dx) {

	PentaDiagonal matrix = setup<PentaDiagonal>(order, dx, coef1::c4);

	// Adjust finite difference approximations at boundary.
	boundary(0, dx, coef1::b1, matrix);
	boundary(1, dx, coef1::b2, matrix);
	boundary(order - 2, dx, coef1::f2, matrix);
	boundary(order - 1, dx, coef1::f1, matrix);

	return matrix;
}

// Forward difference; 1st order accuracy.
TriDiagonal d1dx1::f1(const int order, const double dx) {

	return setup<TriDiagonal>(order, dx, coef1::f1);

	// TODO: Adjust boundary.

}

// Forward difference; 2nd order accuracy.
PentaDiagonal d1dx1::f2(const int order, const double dx) {

	return setup<PentaDiagonal>(order, dx, coef1::f2);

	// TODO: Adjust boundary.

}

// Backward difference; 1st order accuracy.
TriDiagonal d1dx1::b1(const int order, const double dx) {

	return setup<TriDiagonal>(order, dx, coef1::b1);

	// TODO: Adjust boundary.

}

// Backward difference; 2nd order accuracy.
PentaDiagonal d1dx1::b2(const int order, const double dx) {

	return setup<PentaDiagonal>(order, dx, coef1::b2);

	// TODO: Adjust boundary.

}
