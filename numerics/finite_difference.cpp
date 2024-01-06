#include <vector>

#include "band_diagonal_matrix.h"
#include "finite_difference.h"


// Finite difference coefficients.
namespace coefficients {

	// Central difference; 2nd order accuracy.
	std::vector<double> central_2o{ 1.0 / 2.0, 0.0, -1.0 / 2.0 };

	// Central difference; 4th order accuracy.
	std::vector<double> central_4o{ -1.0 / 12.0, 2.0 / 3.0, 0.0, -2.0 / 3.0, 1.0 / 12.0 };

	// Forward difference; 1st order accuracy.
	std::vector<double> forward_1o{ 1.0, -1.0, 0.0 };

	// Forward difference; 2nd order accuracy.
	std::vector<double> forward_2o{ -1.0 / 2.0, 2.0, -3.0 / 2.0, 0.0, 0.0 };

	// Backward difference; 1st order accuracy.
	std::vector<double> backward_1o{ 0.0, 1.0, -1.0 };

	// Backward difference; 2nd order accuracy.
	std::vector<double> backward_2o{ 0.0, 0.0, 3.0 / 2.0, -2.0, 1.0 / 2.0 };

}

// Setting up matrix representation of derivative operator.
template <class T> 
T setup(const int order, const double dx, const std::vector<double>& coef) {

	T matrix(order);

	for (int i = 0; i != coef.size(); ++i)
		for (int j = 0; j != order; ++j)
			matrix.m[i][j] = coef[i] / dx;

	return matrix;

}

// Adjusting finite difference approximations at boundaries.
template <class T>
void boundaries(const int row_index, const double dx, const std::vector<double>& coef, T& matrix) {
	
	for (int i = 0; i != coef.size(); ++i)
			matrix.m[i][row_index] = coef[i] / dx;

}

// Central difference; 2nd order accuracy.
TriDiagonal d1dx1::central_tri(const int order, const double dx) {

	TriDiagonal matrix = setup<TriDiagonal>(order, dx, coefficients::central_2o);

	// Adjust finite difference approximations at boundaries.
	boundaries(0, dx, coefficients::backward_1o, matrix);
	boundaries(order - 1, dx, coefficients::forward_1o, matrix);

	return matrix;

}

// Central difference; 4th order accuracy.
PentaDiagonal d1dx1::central_penta(const int order, const double dx) {

	PentaDiagonal matrix = setup<PentaDiagonal>(order, dx, coefficients::central_4o);

	// Adjust finite difference approximations at boundaries.
	boundaries(0, dx, coefficients::backward_1o, matrix);
	boundaries(1, dx, coefficients::backward_2o, matrix);
	boundaries(order - 2, dx, coefficients::forward_2o, matrix);
	boundaries(order - 1, dx, coefficients::forward_1o, matrix);

	return matrix;
}

// Forward difference; 1st order accuracy.
TriDiagonal d1dx1::forward_1o(const int order, const double dx) {

	return setup<TriDiagonal>(order, dx, coefficients::forward_1o);

	// TODO: Adjust boundaries.

}

// Forward difference; 2nd order accuracy.
PentaDiagonal d1dx1::forward_2o(const int order, const double dx) {

	return setup<PentaDiagonal>(order, dx, coefficients::forward_2o);

	// TODO: Adjust boundaries.

}

// Backward difference; 1st order accuracy.
TriDiagonal d1dx1::backward_1o(const int order, const double dx) {

	return setup<TriDiagonal>(order, dx, coefficients::backward_1o);

	// TODO: Adjust boundaries.

}

// Backward difference; 2nd order accuracy.
PentaDiagonal d1dx1::backward_2o(const int order, const double dx) {

	return setup<PentaDiagonal>(order, dx, coefficients::backward_2o);

	// TODO: Adjust boundaries.

}
