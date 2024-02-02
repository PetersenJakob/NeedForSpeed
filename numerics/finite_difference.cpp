#include <algorithm>
#include <cmath>
#include <stdexcept>
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


// Finite difference coefficients for first order derivative operator.
namespace coef1 {

	// Finite difference representations on equidistant grid.
	// See Fornberg (1988).
	namespace equidistant {

		// Central difference; 2nd order accuracy.
		std::vector<double> c2{
			-1.0 / 2.0,
			 0.0,
			 1.0 / 2.0
		};

		// Central difference; 4th order accuracy.
		std::vector<double> c4{
			 1.0 / 12.0,
			-2.0 / 3.0,
			 0.0,
			 2.0 / 3.0,
			-1.0 / 12.0
		};

		// Forward difference; 1st order accuracy.
		std::vector<double> f1{
			-1.0,
			 1.0
		};

		// Forward difference; 2nd order accuracy.
		std::vector<double> f2{
			-3.0 / 2.0,
			 2.0,
			-1.0 / 2.0
		};

		// Forward difference; 3rd order accuracy.
		std::vector<double> f3{
			-11.0 / 6.0,
			 3.0,
			-3.0 / 2.0,
			 1.0 / 3.0
		};

		// Forward difference; 4th order accuracy.
		std::vector<double> f4{
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

	// Finite difference representations on non-equidistant grid.
	// See Sundqvist and Veronis (1970).
	namespace nonequidistant {

		std::vector<double> c2(
			const double dx_minus,
			const double dx_plus) {

			const double denominator = dx_plus * (1.0 + dx_plus / dx_minus);

			std::vector<double> row(3, 0.0);

			// Sub-diagonal element.
			row[0] = -pow(dx_plus / dx_minus, 2.0) / denominator;

			// Main diagonal element.
			row[1] = -(1.0 - pow(dx_plus / dx_minus, 2.0)) / denominator;

			// Super-diagonal element.
			row[2] = 1.0 / denominator;

			return row;

		}

	}

}


// Finite difference coefficients for second order derivative operator.
namespace coef2 {

	// Finite difference representations on equidistant grid.
	// See Fornberg (1988).
	namespace equidistant {

		// Central difference; 2nd order accuracy.
		std::vector<double> c2{
			 1.0,
			-2.0,
			 1.0
		};

		// Central difference; 4th order accuracy.
		std::vector<double> c4{
			-1.0 / 12.0,
			 4.0 / 3.0,
			-5.0 / 2.0,
			 4.0 / 3.0,
			-1.0 / 12.0
		};

		// Forward difference; 1st order accuracy.
		std::vector<double> f1{
			 1.0,
			-2.0,
			 1.0
		};

		// Forward difference; 2nd order accuracy.
		std::vector<double> f2{
			 2.0,
			-5.0,
			 4.0,
			-1.0
		};

		// Forward difference; 3rd order accuracy.
		std::vector<double> f3{
			 35.0 / 12.0,
			-26.0 / 3.0,
			 19.0 / 2.0,
			-14.0 / 3.0,
			 11.0 / 12.0
		};

		// Forward difference; 4th order accuracy.
		std::vector<double> f4{
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

	// Finite difference representations on non-equidistant grid.
	// See Sundqvist and Veronis (1970).
	namespace nonequidistant {

		std::vector<double> c2(
			const double dx_minus,
			const double dx_plus) {

			const double denominator = dx_plus * dx_minus * (1.0 + dx_plus / dx_minus);

			std::vector<double> row(3, 0.0);

			// Sub-diagonal element.
			row[0] = 2.0 * (dx_plus / dx_minus) / denominator;

			// Main diagonal element.
			row[1] = -2.0 * (1.0 + dx_plus / dx_minus) / denominator;

			// Super-diagonal element.
			row[2] = 2.0 / denominator;

			return row;

		}

	}

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
TriDiagonal d1dx1::equidistant::c2b1(const int order, const double dx) {

	TriDiagonal matrix = setup<TriDiagonal>(order, dx, coef1::equidistant::c2, 1, 2);
	boundary<TriDiagonal>(0, dx, coef1::equidistant::f1, matrix);
	boundary<TriDiagonal>(1, dx, coef1::equidistant::b1, matrix);

	return matrix;

}


// First order derivative operator. 
// Central difference; 2nd order accuracy. Boundary; 2nd order accuracy.
TriDiagonal d1dx1::equidistant::c2b2(const int order, const double dx) {

	TriDiagonal matrix = setup<TriDiagonal>(order, dx, coef1::equidistant::c2, 1, 3);
	boundary<TriDiagonal>(0, dx, coef1::equidistant::f2, matrix);
	boundary<TriDiagonal>(1, dx, coef1::equidistant::b2, matrix);

	return matrix;

}


// First order derivative operator. 
// Central difference; 4th order accuracy. Boundary; 4th order accuracy.
PentaDiagonal d1dx1::equidistant::c4b4(const int order, const double dx) {

	PentaDiagonal matrix = setup<PentaDiagonal>(order, dx, coef1::equidistant::c4, 2, 5);
	boundary<PentaDiagonal>(0, dx, coef1::equidistant::f4, matrix);
	boundary<PentaDiagonal>(2, dx, coef1::equidistant::b4, matrix);
	boundary<PentaDiagonal>(3, dx, coef1::equidistant::b4, matrix);

	return matrix;

}


// Second order derivative operator.
// Central difference; 2nd order accuracy. Boundary; 1st order accuracy.
TriDiagonal d2dx2::equidistant::c2b1(const int order, const double dx) {

	TriDiagonal matrix = setup<TriDiagonal>(order, pow(dx, 2.0), coef2::equidistant::c2, 1, 3);
	boundary<TriDiagonal>(0, pow(dx, 2.0), coef2::equidistant::f1, matrix);
	boundary<TriDiagonal>(1, pow(dx, 2.0), coef2::equidistant::b1, matrix);

	return matrix;

}


// Second order derivative operator.
// Central difference; 2nd order accuracy. Boundary; 2nd order accuracy.
TriDiagonal d2dx2::equidistant::c2b2(const int order, const double dx) {

	TriDiagonal matrix = setup<TriDiagonal>(order, pow(dx, 2.0), coef2::equidistant::c2, 1, 4);
	boundary<TriDiagonal>(0, pow(dx, 2.0), coef2::equidistant::f2, matrix);
	boundary<TriDiagonal>(1, pow(dx, 2.0), coef2::equidistant::b2, matrix);

	return matrix;

}


// Second order derivative operator.
// Central difference; 4th order accuracy. Boundary; 4th order accuracy.
PentaDiagonal d2dx2::equidistant::c4b4(const int order, const double dx) {

	PentaDiagonal matrix = setup<PentaDiagonal>(order, pow(dx, 2.0), coef2::equidistant::c4, 2, 3);
	boundary<PentaDiagonal>(0, pow(dx, 2.0), coef2::equidistant::f1, matrix);
	boundary<PentaDiagonal>(1, pow(dx, 2.0), coef2::equidistant::f1, matrix);
	boundary<PentaDiagonal>(2, pow(dx, 2.0), coef2::equidistant::b1, matrix);
	boundary<PentaDiagonal>(3, pow(dx, 2.0), coef2::equidistant::b1, matrix);

#if false
	PentaDiagonal matrix = setup<PentaDiagonal>(order, pow(dx, 2.0), coef2::c4, 2, 6);
	boundary<PentaDiagonal>(0, pow(dx, 2.0), coef2::f4, matrix);
	boundary<PentaDiagonal>(1, pow(dx, 2.0), coef2::f4, matrix);
	boundary<PentaDiagonal>(2, pow(dx, 2.0), coef2::b4, matrix);
	boundary<PentaDiagonal>(3, pow(dx, 2.0), coef2::b4, matrix);
#endif
	return matrix;

}


// Setting up finite difference representation of derivative operator on non-equidistant grid.
template <class T>
T setup(
	const int order,
	const int derivative_order,
	const std::vector<double> grid,
	const int n_boundary_rows,
	const int n_boundary_elements) {

	T matrix(order, n_boundary_rows, (n_boundary_rows - 1) + n_boundary_elements);

	for (int i = n_boundary_rows; i != order - n_boundary_rows; ++i) {

		double dx_minus = grid[i] - grid[i - 1];
		double dx_plus = grid[i + 1] - grid[i];

		std::vector<double> vec_tmp;
		
		if (std::is_same<T, TriDiagonal>::value) {
			if (derivative_order == 1) {
				vec_tmp = coef1::nonequidistant::c2(dx_minus, dx_plus);
			}
			else if (derivative_order == 2) {
				vec_tmp = coef2::nonequidistant::c2(dx_minus, dx_plus);
			}
			else {
				throw std::invalid_argument("Derivative order should be 1 or 2!");
			}
		}
		else {
			throw std::invalid_argument("Only TriDiagonal for now!!!");
		}

		for (int j = 0; j != matrix.n_diagonals(); ++j) {
			matrix.matrix[j][i] = vec_tmp[j];
		}
	}

	return matrix;

}


// First order derivative operator.
// Central difference; 2nd order accuracy. Boundary; 1st order accuracy.
TriDiagonal d1dx1::nonequidistant::c2b1(const int order, const std::vector<double> grid) {

	TriDiagonal matrix = setup<TriDiagonal>(order, 2, grid, 1, 2);

	const double dx_first = grid[1] - grid[0];
	const double dx_last = grid[order - 1] - grid[order - 2];
	boundary<TriDiagonal>(0, dx_first, coef1::equidistant::f1, matrix);
	boundary<TriDiagonal>(1, dx_last, coef1::equidistant::b1, matrix);

	return matrix;

}
