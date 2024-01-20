#include <algorithm>
#include <cmath>
#include <vector>

#include "band_diagonal_matrix.h"
#include "finite_difference.h"


// Reverse order of coefficients and multiply by scalar.
std::vector<double> reverse_order(std::vector<double> coef, const double scalar = 1.0) {

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
		-1.0 / 2.0,
		 0.0,
		 1.0 / 2.0
	};

    // Central difference; 4th order accuracy.
    std::vector<double> c4 {
		 1.0 / 12.0,
		-2.0 / 3.0,
		 0.0,
		 2.0 / 3.0,
		-1.0 / 12.0
	};

    // Forward difference; 1st order accuracy.
    std::vector<double> f1 {
		-1.0,
		 1.0
	};

    // Forward difference; 2nd order accuracy.
    std::vector<double> f2 {
		-3.0 / 2.0,
		 2.0,
		-1.0 / 2.0
	};

	// Forward difference; 3rd order accuracy.
	std::vector<double> f3 {
		-11.0 / 6.0,
		 3.0,
		-3.0 / 2.0,
		 1.0 / 3.0
	};

	// Forward difference; 4th order accuracy.
	std::vector<double> f4 {
		-25.0 / 12.0,
		 4.0,
		-3.0,
		 4.0 / 3.0,
		-1.0 / 4.0
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
		 1.0
	};

	// Forward difference; 2nd order accuracy.
	std::vector<double> f2 {
		 2.0,
		-5.0,
		 4.0,
		-1.0
	};

	// Forward difference; 3rd order accuracy.
	std::vector<double> f3 {
		 35.0 / 12.0,
		-26.0 / 3.0,
		 19.0 / 2.0,
		-14.0 / 3.0,
		 11.0 / 12.0
	};

	// Forward difference; 4th order accuracy.
	std::vector<double> f4 {
		 15.0 / 4.0,
		-77.0 / 6.0,
		 107.0 / 6.0,
		-13.0,
		 61.0 / 12.0,
		-5.0 / 6.0
	};

	// Backward difference; 1st order accuracy.
	std::vector<double> b1 = reverse_order(f1);

	// Backward difference; 2nd order accuracy.
	std::vector<double> b2 = reverse_order(f2);

	// Backward difference; 3rd order accuracy.
	std::vector<double> b3 = reverse_order(f3);

	// Backward difference; 4th order accuracy.
	std::vector<double> b4 = reverse_order(f4);

}


// Setting up finite difference representation of derivative operator on equidistant grid.
template <class T> 
T setup(
	const int order, 
	const double dx, 
	const std::vector<double>& coef,
	const int n_boundary_rows,
	const int n_boundary_elements) {

	T matrix(order, n_boundary_rows, (n_boundary_rows - 1) + n_boundary_elements);

	for (int i = 0; i != coef.size(); ++i) {
		for (int j = 0; j != order; ++j) {
			matrix.matrix[i][j] = coef[i] / dx;
		}
	}

	return matrix;

}


// Adjusting finite difference representation at boundary.
template <class T>
void boundary(const int row_index, const double dx, const std::vector<double>& coef, T& matrix) {
	
	int index_tmp = row_index - matrix.n_boundary_rows();
	if (index_tmp < 0) {
		index_tmp = row_index;
	}

	for (int i = 0; i != coef.size(); ++i) {
		matrix.boundary_rows[row_index][index_tmp + i] = coef[i] / dx;
	}

}


// First order derivative operator. 
// Central difference; 2nd order accuracy. Boundary; 1st order accuracy.
TriDiagonal d1dx1::c2b1(const int order, const double dx) {

	TriDiagonal matrix = setup<TriDiagonal>(order, dx, coef1::c2, 1, 2);
	boundary(0, dx, coef1::f1, matrix);
	boundary(1, dx, coef1::b1, matrix);

	return matrix;

}


// First order derivative operator. 
// Central difference; 2nd order accuracy. Boundary; 2nd order accuracy.
TriDiagonal d1dx1::c2b2(const int order, const double dx) {

	TriDiagonal matrix = setup<TriDiagonal>(order, dx, coef1::c2, 1, 3);
	boundary(0, dx, coef1::f2, matrix);
	boundary(1, dx, coef1::b2, matrix);

	return matrix;

}


// First order derivative operator. 
// Central difference; 4th order accuracy. Boundary; 4th order accuracy.
PentaDiagonal d1dx1::c4b4(const int order, const double dx) {

	PentaDiagonal matrix = setup<PentaDiagonal>(order, dx, coef1::c4, 2, 5);
	boundary(0, dx, coef1::f4, matrix);
	boundary(1, dx, coef1::f4, matrix);
	boundary(2, dx, coef1::b4, matrix);
	boundary(3, dx, coef1::b4, matrix);

	return matrix;

}


// Second order derivative operator.
// Central difference; 2nd order accuracy. Boundary; 1st order accuracy.
TriDiagonal d2dx2::c2b1(const int order, const double dx) {

	TriDiagonal matrix = setup<TriDiagonal>(order, pow(dx, 2.0), coef2::c2, 1, 3);
	boundary(0, pow(dx, 2.0), coef2::f1, matrix);
	boundary(1, pow(dx, 2.0), coef2::b1, matrix);

	std::vector<double> coefficients(3, 0.0);
//	boundary(0, pow(dx, 2.0), coefficients, matrix);
//	boundary(1, pow(dx, 2.0), coefficients, matrix);

	return matrix;

}


// Second order derivative operator.
// Central difference; 2nd order accuracy. Boundary; 2nd order accuracy.
TriDiagonal d2dx2::c2b2(const int order, const double dx) {

	TriDiagonal matrix = setup<TriDiagonal>(order, pow(dx, 2.0), coef2::c2, 1, 4);
	boundary(0, pow(dx, 2.0), coef2::f2, matrix);
	boundary(1, pow(dx, 2.0), coef2::b2, matrix);

	return matrix;

}


// Second order derivative operator.
// Central difference; 4th order accuracy. Boundary; 4th order accuracy.
PentaDiagonal d2dx2::c4b4(const int order, const double dx) {

	PentaDiagonal matrix = setup<PentaDiagonal>(order, pow(dx, 2.0), coef2::c4, 2, 6);
	boundary(0, pow(dx, 2.0), coef2::f4, matrix);
	boundary(1, pow(dx, 2.0), coef2::f4, matrix);
	boundary(2, pow(dx, 2.0), coef2::b4, matrix);
	boundary(3, pow(dx, 2.0), coef2::b4, matrix);

	return matrix;

}
