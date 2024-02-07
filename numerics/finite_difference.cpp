#include <algorithm>
#include <cmath>
#include <functional>
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


// ...
std::vector<double> adjust_coefficients(std::vector<double> coefficients, const double denominator) {

	for (auto& element : coefficients) {
		element /= denominator;
	}

	return coefficients;

}


// Finite difference coefficients for first order derivative operator.
namespace coef1 {

	// Finite difference representations on equidistant grid.
	// See Fornberg (1988).
	namespace equidistant {

		// Central difference; 2nd order accuracy.
		std::vector<double> c2_coefficients{

			-1.0 / 2.0,
			 0.0,
			 1.0 / 2.0

		};

		std::vector<double> c2(const double dx) {

			return adjust_coefficients(c2_coefficients, dx);

		}

		// Central difference; 4th order accuracy.
		std::vector<double> c4_coefficients{

			 1.0 / 12.0,
			-2.0 / 3.0,
			 0.0,
			 2.0 / 3.0,
			-1.0 / 12.0

		};

		std::vector<double> c4(const double dx) {

			return adjust_coefficients(c4_coefficients, dx);

		}

		// Forward difference; 1st order accuracy.
		std::vector<double> f1_coefficients{

			-1.0,
			 1.0

		};

		std::vector<double> f1(const double dx) {

			return adjust_coefficients(f1_coefficients, dx);

		}

		// Forward difference; 2nd order accuracy.
		std::vector<double> f2_coefficients{

			-3.0 / 2.0,
			 2.0,
			-1.0 / 2.0

		};

		std::vector<double> f2(const double dx) {

			return adjust_coefficients(f2_coefficients, dx);

		}

		// Forward difference; 3rd order accuracy.
		std::vector<double> f3_coefficients{

			-11.0 / 6.0,
			 3.0,
			-3.0 / 2.0,
			 1.0 / 3.0

		};

		std::vector<double> f3(const double dx) {

			return adjust_coefficients(f3_coefficients, dx);

		}

		// Forward difference; 4th order accuracy.
		std::vector<double> f4_coefficients{

			-25.0 / 12.0,
			 4.0,
			-3.0,
			 4.0 / 3.0,
			-1.0 / 4.0

		};

		std::vector<double> f4(const double dx) {

			return adjust_coefficients(f4_coefficients, dx);

		}

		// Backward difference; 1st order accuracy.
		std::vector<double> b1(const double dx) {

			return reverse_order(f1(dx), -1.0);

		}

		// Backward difference; 2nd order accuracy.
		std::vector<double> b2(const double dx) {

			return reverse_order(f2(dx), -1.0);

		}

		// Backward difference; 3rd order accuracy.
		std::vector<double> b3(const double dx) {

			return reverse_order(f3(dx), -1.0);

		}

		// Backward difference; 4th order accuracy.
		std::vector<double> b4(const double dx) {

			return reverse_order(f4(dx), -1.0);

		}

	}

	// Finite difference representations on non-equidistant grid.
	// See Sundqvist and Veronis (1970).
	namespace nonequidistant {

		// dx_vector: Step size vector with four elements [dx(-2), dx(-1), dx(+1), dx(+2)].

		// Central difference; 2nd order accuracy.
		std::vector<double> c2(const std::vector<double>& dx_vector) {

			const double dx_m = dx_vector[1];
			const double dx_p = dx_vector[2];

			std::vector<double> row(3, 0.0);

			const double denominator = dx_p * (1.0 + dx_p / dx_m);

			// Coefficient of 1st sub-diagonal.
			row[0] = -pow(dx_p / dx_m, 2) / denominator;

			// Coefficient of main diagonal.
			row[1] = -(1.0 - pow(dx_p / dx_m, 2)) / denominator;

			// Coefficient of 1st super-diagonal.
			row[2] = 1.0 / denominator;

			return row;

		}

		// Central difference; 4th order accuracy.
		std::vector<double> c4(const std::vector<double>& dx_vector) {

			const double dx_m2 = dx_vector[0];
			const double dx_m1 = dx_vector[1];
			const double dx_p1 = dx_vector[2];
			const double dx_p2 = dx_vector[3];

			std::vector<double> row(5, 0.0);

			const double denominator =
				  pow(dx_m1 + dx_m2, 2) * (dx_p1 + dx_p2)
				- 32.0 * dx_p1 * pow(dx_m1, 2)
				- 32.0 * pow(dx_p1, 2) * dx_m1
				+ (dx_m1 + dx_m2) * pow(dx_p1 + dx_p2, 2);

			// Coefficient of 2nd sub-diagonal.
			row[0] = -pow(dx_p1 + dx_p2, 2) / denominator;

			// Coefficient of 1st sub-diagonal.
			row[1] = 32.0 * pow(dx_p1, 2) / denominator;

			// Coefficient of main diagonal.
			row[2] = -(pow(dx_m1 + dx_m2, 2) - 32.0 * pow(dx_m1, 2) + 32.0 * pow(dx_p1, 2) - pow(dx_p1 + dx_p2, 2)) / denominator;

			// Coefficient of 1st super-diagonal.
			row[3] = -32.0 * pow(dx_m1, 2) / denominator;

			// Coefficient of 2nd super-diagonal.
			row[4] = pow(dx_m1 + dx_m2, 2) / denominator;

			return row;

		}

		// Forward difference; 1st order accuracy.
		std::vector<double> f1(const std::vector<double>& dx_vector) {

			const double dx_p = dx_vector[2];

			std::vector<double> row(2, 0.0);

			const double denominator = dx_p;

			// Coefficient of main diagonal.
			row[0] = -1.0 / denominator;

			// Coefficient of 1st super-diagonal.
			row[1] = 1.0 / denominator;

			return row;

		}

		// Forward difference; 2nd order accuracy.
		std::vector<double> f2(const std::vector<double>& dx_vector) {

			const double dx_p1 = dx_vector[2];
			const double dx_p2 = dx_vector[3];

			std::vector<double> row(3, 0.0);

			const double denominator = pow(dx_p1, 2) * (dx_p1 + dx_p2) - dx_p1 * pow(dx_p1 + dx_p2, 2);;

			// Coefficient of main diagonal.
			row[0] = (pow(dx_p2, 2) + 2 * dx_p1 * dx_p2) / denominator;

			// Coefficient of 1st super-diagonal.
			row[1] = -pow(dx_p1 + dx_p2, 2) / denominator;

			// Coefficient of 2nd super-diagonal.
			row[2] = pow(dx_p1, 2) / denominator;

			return row;

		}

		// TODO: For backward differenc, reverse order of dx_vector??
		// Backward difference; 1st order accuracy.
		std::vector<double> b1(const std::vector<double>& dx_vector) {

			return reverse_order(f1(dx_vector), -1.0);

		}

		// Backward difference; 2nd order accuracy.
		std::vector<double> b2(const std::vector<double>& dx_vector) {

			return reverse_order(f2(dx_vector), -1.0);

		}

	}

}


// Finite difference coefficients for second order derivative operator.
namespace coef2 {

	// Finite difference representations on equidistant grid.
	// See Fornberg (1988).
	namespace equidistant {

		// Central difference; 2nd order accuracy.
		std::vector<double> c2_coefficients{

			 1.0,
			-2.0,
			 1.0

		};

		std::vector<double> c2(const double dx) {

			return adjust_coefficients(c2_coefficients, dx * dx);

		}

		// Central difference; 4th order accuracy.
		std::vector<double> c4_coefficients{

			-1.0 / 12.0,
			 4.0 / 3.0,
			-5.0 / 2.0,
			 4.0 / 3.0,
			-1.0 / 12.0

		};

		std::vector<double> c4(const double dx) {

			return adjust_coefficients(c4_coefficients, dx * dx);

		}

		// Forward difference; 1st order accuracy.
		std::vector<double> f1_coefficients{

			 1.0,
			-2.0,
			 1.0

		};

		std::vector<double> f1(const double dx) {

			return adjust_coefficients(f1_coefficients, dx * dx);

		}

		// Forward difference; 2nd order accuracy.
		std::vector<double> f2_coefficients{

			 2.0,
			-5.0,
			 4.0,
			-1.0

		};
		
		std::vector<double> f2(const double dx) {

			return adjust_coefficients(f2_coefficients, dx * dx);

		}

		// Forward difference; 3rd order accuracy.
		std::vector<double> f3_coefficients{

			 35.0 / 12.0,
			-26.0 / 3.0,
			 19.0 / 2.0,
			-14.0 / 3.0,
			 11.0 / 12.0

		};

		std::vector<double> f3(const double dx) {

			return adjust_coefficients(f3_coefficients, dx * dx);

		}

		// Forward difference; 4th order accuracy.
		std::vector<double> f4_coefficients{

			 15.0 / 4.0,
			-77.0 / 6.0,
			 107.0 / 6.0,
			-13.0,
			 61.0 / 12.0,
			-5.0 / 6.0

		};

		std::vector<double> f4(const double dx) {

			return adjust_coefficients(f4_coefficients, dx * dx);

		}

		// Backward difference; 1st order accuracy.
		std::vector<double> b1(const double dx) {

			return reverse_order(f1(dx));

		}

		// Backward difference; 2nd order accuracy.
		std::vector<double> b2(const double dx) {

			return reverse_order(f2(dx));

		}

		// Backward difference; 3rd order accuracy.
		std::vector<double> b3(const double dx) {

			return reverse_order(f3(dx));

		}

		// Backward difference; 4th order accuracy.
		std::vector<double> b4(const double dx) {

			return reverse_order(f4(dx));

		}

	}

	// Finite difference representations on non-equidistant grid.
	// See Sundqvist and Veronis (1970).
	namespace nonequidistant {

		// dx_vector: Step size vector with four elements [dx(-2), dx(-1), dx(+1), dx(+2)].

		// Central difference; ~2nd order accuracy.
		std::vector<double> c2(const std::vector<double>& dx_vector) {

			const double dx_m = dx_vector[1];
			const double dx_p = dx_vector[2];

			std::vector<double> row(3, 0.0);

			const double denominator = dx_p * dx_m * (1.0 + dx_p / dx_m);

			// Sub-diagonal element.
			row[0] = 2.0 * (dx_p / dx_m) / denominator;

			// Main diagonal element.
			row[1] = -2.0 * (1.0 + dx_p / dx_m) / denominator;

			// Super-diagonal element.
			row[2] = 2.0 / denominator;

			return row;

		}

		// Central difference; ~4th order accuracy.
		std::vector<double> c4(const std::vector<double>& dx_vector) {

			const double dx_m2 = dx_vector[0];
			const double dx_m1 = dx_vector[1];
			const double dx_p1 = dx_vector[2];
			const double dx_p2 = dx_vector[3];

			std::vector<double> row(5, 0.0);

			const double denominator =
				- (dx_m1 + dx_m2) * pow(dx_p1 + dx_p2, 2) / 2.0
				+ 16.0 * pow(dx_p1, 2) * dx_m1
				+ 16.0 * dx_p1 * pow(dx_m1, 2)
				- pow(dx_m1 + dx_m2, 2) * (dx_p1 + dx_p2) / 2.0;

			// Coefficient of 2nd sub-diagonal.
			row[0] = -(dx_p1 + dx_p2) / denominator;

			// Coefficient of 1st sub-diagonal.
			row[1] = 32.0 * dx_p1 / denominator;

			// Coefficient of main diagonal.
			row[2] = -(-(dx_m1 + dx_m2) + 32.0 * dx_m1 + 32.0 * dx_p1 - (dx_p1 + dx_p2)) / denominator;

			// Coefficient of 1st super-diagonal.
			row[3] = 32.0 * dx_m1 / denominator;

			// Coefficient of 2nd super-diagonal.
			row[4] = -(dx_m1 + dx_m2) / denominator;

			return row;

		}

	}

}


// Setting up finite difference representation of derivative operator on equidistant grid.
template <class T>
T setup(
	const int order,
	const std::vector<double>& coef,
	const int n_boundary_rows,
	const int n_boundary_elements) {

	T matrix(order, n_boundary_rows, (n_boundary_rows - 1) + n_boundary_elements);

	for (int i = 0; i != coef.size(); ++i) {

		for (int j = n_boundary_rows; j != order - n_boundary_rows; ++j) {
//		for (int j = 0; j != order; ++j) {
			matrix.matrix[i][j] = coef[i];
		}
	}

	return matrix;

}


// Adjusting finite difference representation at boundary.
template <class T>
void boundary(const int row_index, const std::vector<double>& coef, T& matrix) {

	int index_tmp = row_index - matrix.n_boundary_rows();
	if (index_tmp < 0) {
		index_tmp = row_index;
	}

	for (int i = 0; i != coef.size(); ++i) {
		matrix.boundary_rows[row_index][index_tmp + i] = coef[i];
	}

}


// First order derivative operator. 
// Central difference; 2nd order accuracy. Boundary; 1st order accuracy.
TriDiagonal d1dx1::equidistant::c2b1(const int order, const double dx) {

	// TODO: If you choose two boundary rows with f1+f1 and b1+b1, what will the L2 function norm be?

	TriDiagonal matrix = setup<TriDiagonal>(order, coef1::equidistant::c2(dx), 1, 2);
	boundary<TriDiagonal>(0, coef1::equidistant::f1(dx), matrix);
	boundary<TriDiagonal>(1, coef1::equidistant::b1(dx), matrix);

	return matrix;

}


// First order derivative operator. 
// Central difference; 2nd order accuracy. Boundary; 2nd order accuracy.
TriDiagonal d1dx1::equidistant::c2b2(const int order, const double dx) {

	TriDiagonal matrix = setup<TriDiagonal>(order, coef1::equidistant::c2(dx), 1, 3);
	boundary<TriDiagonal>(0, coef1::equidistant::f2(dx), matrix);
	boundary<TriDiagonal>(1, coef1::equidistant::b2(dx), matrix);

	return matrix;

}


// First order derivative operator. 
// Central difference; 4th order accuracy. Boundary; 2nd order accuracy.
PentaDiagonal d1dx1::equidistant::c4b2(const int order, const double dx) {

	PentaDiagonal matrix = setup<PentaDiagonal>(order, coef1::equidistant::c4(dx), 2, 3);
	boundary<PentaDiagonal>(0, coef1::equidistant::f2(dx), matrix);
	boundary<PentaDiagonal>(1, coef1::equidistant::f2(dx), matrix);
	boundary<PentaDiagonal>(2, coef1::equidistant::b2(dx), matrix);
	boundary<PentaDiagonal>(3, coef1::equidistant::b2(dx), matrix);

	return matrix;

}


// First order derivative operator. 
// Central difference; 4th order accuracy. Boundary; 4th order accuracy.
PentaDiagonal d1dx1::equidistant::c4b4(const int order, const double dx) {

	PentaDiagonal matrix = setup<PentaDiagonal>(order, coef1::equidistant::c4(dx), 2, 5);
	boundary<PentaDiagonal>(0, coef1::equidistant::f4(dx), matrix);
	boundary<PentaDiagonal>(1, coef1::equidistant::f4(dx), matrix);
	boundary<PentaDiagonal>(2, coef1::equidistant::b4(dx), matrix);
	boundary<PentaDiagonal>(3, coef1::equidistant::b4(dx), matrix);

	return matrix;

}


// Second order derivative operator.
// Central difference; 2nd order accuracy. Boundary; 1st order accuracy.
TriDiagonal d2dx2::equidistant::c2b1(const int order, const double dx) {

	TriDiagonal matrix = setup<TriDiagonal>(order, coef2::equidistant::c2(dx), 1, 3);
	boundary<TriDiagonal>(0, coef2::equidistant::f1(dx), matrix);
	boundary<TriDiagonal>(1, coef2::equidistant::b1(dx), matrix);

	return matrix;

}


// Second order derivative operator.
// Central difference; 2nd order accuracy. Boundary; 2nd order accuracy.
TriDiagonal d2dx2::equidistant::c2b2(const int order, const double dx) {

	TriDiagonal matrix = setup<TriDiagonal>(order, coef2::equidistant::c2(dx), 1, 4);
	boundary<TriDiagonal>(0, coef2::equidistant::f2(dx), matrix);
	boundary<TriDiagonal>(1, coef2::equidistant::b2(dx), matrix);

	return matrix;

}


// Second order derivative operator.
// Central difference; 4th order accuracy. Boundary; TODO 1st order accuracy.
PentaDiagonal d2dx2::equidistant::c4b4(const int order, const double dx) {

#if false
	PentaDiagonal matrix = setup<PentaDiagonal>(order, coef2::equidistant::c4(dx), 2, 3);
	boundary<PentaDiagonal>(0, coef2::equidistant::f1(dx), matrix);
	boundary<PentaDiagonal>(1, coef2::equidistant::f1(dx), matrix);
	boundary<PentaDiagonal>(2, coef2::equidistant::b1(dx), matrix);
	boundary<PentaDiagonal>(3, coef2::equidistant::b1(dx), matrix);
#endif

	PentaDiagonal matrix = setup<PentaDiagonal>(order, coef2::equidistant::c4(dx), 2, 6);
	boundary<PentaDiagonal>(0, coef2::equidistant::f4(dx), matrix);
	boundary<PentaDiagonal>(1, coef2::equidistant::f4(dx), matrix);
	boundary<PentaDiagonal>(2, coef2::equidistant::b4(dx), matrix);
	boundary<PentaDiagonal>(3, coef2::equidistant::b4(dx), matrix);

	return matrix;

}


// Setting up finite difference representation of derivative operator on non-equidistant grid.
template <class T>
T setup(
	const int order,
	const std::vector<double> grid,
	std::function<std::vector<double>(std::vector<double>)> coef,
	const int n_boundary_rows,
	const int n_boundary_elements) {

	T matrix(order, n_boundary_rows, (n_boundary_rows - 1) + n_boundary_elements);

	std::vector<double> coef_tmp;

	for (int i = n_boundary_rows; i != order - n_boundary_rows; ++i) {

		if (std::is_same<T, TriDiagonal>::value) {

			int idx_m1 = i - 1;
			int idx_p1 = i + 1;

			double dx_m1 = grid[i] - grid[idx_m1];
			double dx_p1 = grid[idx_p1] - grid[i];

			coef_tmp = coef({ 0.0, dx_m1, dx_p1, 0.0 });

		}
		else if (std::is_same<T, PentaDiagonal>::value) {

			int idx_m2 = i - 2;
			int idx_m1 = i - 1;
			int idx_p1 = i + 1;
			int idx_p2 = i + 2;

			double dx_m2 = grid[idx_m1] - grid[idx_m2];
			double dx_m1 = grid[i] - grid[idx_m1];
			double dx_p1 = grid[idx_p1] - grid[i];
			double dx_p2 = grid[idx_p2] - grid[idx_p1];

			coef_tmp = coef({ dx_m2, dx_m1, dx_p1, dx_p2 });

		}
		else {
			throw std::invalid_argument("Unknown matrix.");
		}

		for (int j = 0; j != matrix.n_diagonals(); ++j) {
			matrix.matrix[j][i] = coef_tmp[j];
		}

	}

	return matrix;

}


// First order derivative operator.
// Central difference; 2nd order accuracy. Boundary; 1st order accuracy.
TriDiagonal d1dx1::nonequidistant::c2b1(const int order, const std::vector<double> grid) {

	TriDiagonal matrix = setup<TriDiagonal>(order, grid, coef1::nonequidistant::c2, 1, 2);

	std::vector<double> dx_vec_first = { 0.0, 0.0, grid[1] - grid[0], 0.0 };
	boundary<TriDiagonal>(0, coef1::nonequidistant::f1(dx_vec_first), matrix);
	// TODO: Should dx_last be reversed?
	std::vector<double> dx_vec_last = { 0.0, 0.0, grid[order - 1] - grid[order - 2], 0.0 };
	boundary<TriDiagonal>(1, coef1::nonequidistant::b1(dx_vec_last), matrix);

	return matrix;

}


// First order derivative operator.
// Central difference; 2nd order accuracy. Boundary; 2nd order accuracy.
TriDiagonal d1dx1::nonequidistant::c2b2(const int order, const std::vector<double> grid) {

	TriDiagonal matrix = setup<TriDiagonal>(order, grid, coef1::nonequidistant::c2, 1, 3);

	std::vector<double> dx_vec_first = { 0.0, 0.0, grid[1] - grid[0], grid[2] - grid[1] };
	boundary<TriDiagonal>(0, coef1::nonequidistant::f2(dx_vec_first), matrix);

	// TODO: Should dx_last be reversed? Yes!
	std::vector<double> dx_vec_last = { 0.0, 0.0, grid[order - 1] - grid[order - 2], grid[order - 2] - grid[order - 3] };
	boundary<TriDiagonal>(1, coef1::nonequidistant::b2(dx_vec_last), matrix);

	return matrix;

}


// First order derivative operator.
// Central difference; 4th order accuracy. Boundary; 2nd order accuracy.
PentaDiagonal d1dx1::nonequidistant::c4b2(const int order, const std::vector<double> grid) {

	// TODO: Be able to choose c2 as second boundary row, and f1 as first boundary row.

	PentaDiagonal matrix = setup<PentaDiagonal>(order, grid, coef1::nonequidistant::c4, 2, 3);

	std::vector<double> dx_vec_1 = { 0.0, 0.0, grid[1] - grid[0], grid[2] - grid[1] };
	boundary<PentaDiagonal>(0, coef1::nonequidistant::f2(dx_vec_1), matrix);

	std::vector<double> dx_vec_2 = { 0.0, 0.0, grid[2] - grid[1], grid[3] - grid[2] };
	boundary<PentaDiagonal>(1, coef1::nonequidistant::f2(dx_vec_2), matrix);

	std::vector<double> dx_vec_3 = { 0.0, 0.0, grid[order - 2] - grid[order - 3], grid[order - 3] - grid[order - 4] };
	boundary<PentaDiagonal>(2, coef1::nonequidistant::b2(dx_vec_3), matrix);

	std::vector<double> dx_vec_4 = { 0.0, 0.0, grid[order - 1] - grid[order - 2], grid[order - 2] - grid[order - 3] };
	boundary<PentaDiagonal>(3, coef1::nonequidistant::b2(dx_vec_4), matrix);

	return matrix;

}


// Second order derivative operator.
// Central difference; 2nd order accuracy. Boundary; 1st order accuracy.
TriDiagonal d2dx2::nonequidistant::c2b1(const int order, const std::vector<double> grid) {

	TriDiagonal matrix = setup<TriDiagonal>(order, grid, coef2::nonequidistant::c2, 1, 2);

	// f1 and b1 for coef2::nonequidistant...

	std::vector<double> dx_vec_first = { 0.0, 0.0, grid[1] - grid[0], 0.0 };
	boundary<TriDiagonal>(0, coef1::nonequidistant::f1(dx_vec_first), matrix);
	// TODO: Should dx_last be reversed?
	std::vector<double> dx_vec_last = { 0.0, 0.0, grid[order - 1] - grid[order - 2], 0.0 };
	boundary<TriDiagonal>(1, coef1::nonequidistant::b1(dx_vec_last), matrix);

	return matrix;

}
