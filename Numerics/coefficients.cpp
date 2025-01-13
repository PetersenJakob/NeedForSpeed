#include <algorithm>
#include <vector>

#include "coefficients.h"

#if false
// Finite difference coefficients for first order derivative operator.
namespace coef_x1Template {

	// Finite difference representation on uniform grid.
	// See Fornberg (1988).
	namespace uniform {

		// Central difference; 2nd order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> c2_coefficients{
			-1.0 / 2.0,
			 0.0,
			 1.0 / 2.0
		};

		template<typename Tnumber>
		const std::vector<Tnumber> c2(const Tnumber dx) {

			return adjust_coefficients(c2_coefficients<Tnumber>, dx);

		}

		// Central difference; 4th order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> c4_coefficients{
			 1.0 / 12.0,
			-2.0 / 3.0,
			 0.0,
			 2.0 / 3.0,
			-1.0 / 12.0
		};

		template<typename Tnumber>
		const std::vector<Tnumber> c4(const Tnumber dx) {

			return adjust_coefficients(c4_coefficients<Tnumber>, dx);

		}

		// Forward difference; 1st order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> f1_coefficients{
			-1.0,
			 1.0
		};

		template<typename Tnumber>
		const std::vector<Tnumber> f1(const Tnumber dx) {

			return adjust_coefficients(f1_coefficients<Tnumber>, dx);

		}

		// Forward difference; 2nd order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> f2_coefficients{
			-3.0 / 2.0,
			 2.0,
			-1.0 / 2.0
		};

		template<typename Tnumber>
		const std::vector<Tnumber> f2(const Tnumber dx) {

			return adjust_coefficients(f2_coefficients<Tnumber>, dx);

		}

		// Forward difference; 3rd order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> f3_coefficients{
			-11.0 / 6.0,
			 3.0,
			-3.0 / 2.0,
			 1.0 / 3.0
		};

		template<typename Tnumber>
		const std::vector<Tnumber> f3(const Tnumber dx) {

			return adjust_coefficients(f3_coefficients<Tnumber>, dx);

		}

		// Forward difference; 4th order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> f4_coefficients{
			-25.0 / 12.0,
			 4.0,
			-3.0,
			 4.0 / 3.0,
			-1.0 / 4.0
		};

		template<typename Tnumber>
		const std::vector<Tnumber> f4(const Tnumber dx) {

			return adjust_coefficients(f4_coefficients<Tnumber>, dx);

		}

		// Backward difference; 1st order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> b1(const Tnumber dx) {

			return reverse_order(f1(dx), -1.0);

		}

		// Backward difference; 2nd order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> b2(const Tnumber dx) {

			return reverse_order(f2(dx), -1.0);

		}

		// Backward difference; 3rd order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> b3(const Tnumber dx) {

			return reverse_order(f3(dx), -1.0);

		}

		// Backward difference; 4th order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> b4(const Tnumber dx) {

			return reverse_order(f4(dx), -1.0);

		}

	}

	// Finite difference representation on non-uniform grid.
	// See Sundqvist and Veronis (1970).
	namespace nonuniform {

		// dx_vector; step size vector with four elements:
		// [0] dx_m2 = x_m1 - x_m2  TODO: Correct sign?
		// [1] dx_m1 = x - x_m1     TODO: Correct sign?
		// [2] dx_p1 = x_p1 - x  
		// [3] dx_p2 = x_p2 - x_p1

		// Central difference; 2nd order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> c2(const std::vector<Tnumber>& dx_vector) {

			const Tnumber dx_m1 = dx_vector[1];
			const Tnumber dx_p1 = dx_vector[2];

			std::vector<Tnumber> row(3, 0.0);

			const Tnumber denominator = dx_p1 * (1.0 + dx_p1 / dx_m1);

			// Coefficient of 1st sub-diagonal.
			row[0] = -pow(dx_p1 / dx_m1, 2) / denominator;

			// Coefficient of main diagonal.
			row[1] = -(1.0 - pow(dx_p1 / dx_m1, 2)) / denominator;

			// Coefficient of 1st super-diagonal.
			row[2] = 1.0 / denominator;

			return row;

		}

		// Central difference; 4th order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> c4(const std::vector<Tnumber>& dx_vector) {

			const Tnumber dx_m2 = dx_vector[0];
			const Tnumber dx_m1 = dx_vector[1];
			const Tnumber dx_p1 = dx_vector[2];
			const Tnumber dx_p2 = dx_vector[3];

			std::vector<Tnumber> row(5, 0.0);

			const Tnumber denominator =
				pow(dx_m1 + dx_m2, 2) * (dx_p1 + dx_p2)
				- 32.0 * dx_p1 * pow(dx_m1, 2)
				- 32.0 * pow(dx_p1, 2) * dx_m1
				+ (dx_m1 + dx_m2) * pow(dx_p1 + dx_p2, 2);

			// Coefficient of 2nd sub-diagonal.
			row[0] = -pow(dx_p1 + dx_p2, 2) / denominator;

			// Coefficient of 1st sub-diagonal.
			row[1] = 32.0 * pow(dx_p1, 2) / denominator;

			// Coefficient of main diagonal.
			row[2] =
				-(pow(dx_m1 + dx_m2, 2) - 32.0 * pow(dx_m1, 2)
					+ 32.0 * pow(dx_p1, 2) - pow(dx_p1 + dx_p2, 2)) / denominator;

			// Coefficient of 1st super-diagonal.
			row[3] = -32.0 * pow(dx_m1, 2) / denominator;

			// Coefficient of 2nd super-diagonal.
			row[4] = pow(dx_m1 + dx_m2, 2) / denominator;

			return row;

		}

		// Forward difference; 1st order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> f1(const std::vector<Tnumber>& dx_vector) {

			const Tnumber dx_p1 = dx_vector[2];

			std::vector<Tnumber> row(2, 0.0);

			const Tnumber denominator = dx_p1;

			// Coefficient of main diagonal.
			row[0] = -1.0 / denominator;

			// Coefficient of 1st super-diagonal.
			row[1] = 1.0 / denominator;

			return row;

		}

		// Forward difference; 2nd order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> f2(const std::vector<Tnumber>& dx_vector) {

			const Tnumber dx_p1 = dx_vector[2];
			const Tnumber dx_p2 = dx_vector[3];

			std::vector<Tnumber> row(3, 0.0);

			const Tnumber denominator =
				pow(dx_p1, 2) * (dx_p1 + dx_p2) - dx_p1 * pow(dx_p1 + dx_p2, 2);

			// Coefficient of main diagonal.
			row[0] = (pow(dx_p2, 2) + 2 * dx_p1 * dx_p2) / denominator;

			// Coefficient of 1st super-diagonal.
			row[1] = -pow(dx_p1 + dx_p2, 2) / denominator;

			// Coefficient of 2nd super-diagonal.
			row[2] = pow(dx_p1, 2) / denominator;

			return row;

		}

		// Backward difference; 1st order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> b1(const std::vector<Tnumber>& dx_vector) {

			return reverse_order(f1(dx_vector), -1.0);

		}

		// Backward difference; 2nd order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> b2(const std::vector<Tnumber>& dx_vector) {

			return reverse_order(f2(dx_vector), -1.0);

		}

	}

}
#endif

// Finite difference coefficients for second order derivative operator.
namespace coef_x2Template {
	// Finite difference representation on uniform grid.
	// See Fornberg (1988).
	namespace uniform {

		// Central difference; 2nd order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> c2_coefficients{
			 1.0,
			-2.0,
			 1.0
		};

		template<typename Tnumber>
		const std::vector<Tnumber> c2(const Tnumber dx) {

			return adjust_coefficients(c2_coefficients<Tnumber>, dx * dx);

		}

		// Central difference; 4th order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> c4_coefficients{
			-1.0 / 12.0,
			 4.0 / 3.0,
			-5.0 / 2.0,
			 4.0 / 3.0,
			-1.0 / 12.0
		};

		template<typename Tnumber>
		const std::vector<Tnumber> c4(const Tnumber dx) {

			return adjust_coefficients(c4_coefficients<Tnumber>, dx * dx);

		}

		// Forward difference; 1st order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> f1_coefficients{
			 1.0,
			-2.0,
			 1.0
		};

		template<typename Tnumber>
		const std::vector<Tnumber> f1(const Tnumber dx) {

			return adjust_coefficients(f1_coefficients<Tnumber>, dx * dx);

		}

		// Forward difference; 2nd order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> f2_coefficients{
			 2.0,
			-5.0,
			 4.0,
			-1.0
		};

		template<typename Tnumber>
		const std::vector<Tnumber> f2(const Tnumber dx) {

			return adjust_coefficients(f2_coefficients<Tnumber>, dx * dx);

		}

		// Forward difference; 3rd order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> f3_coefficients{
			 35.0 / 12.0,
			-26.0 / 3.0,
			 19.0 / 2.0,
			-14.0 / 3.0,
			 11.0 / 12.0
		};

		template<typename Tnumber>
		const std::vector<Tnumber> f3(const Tnumber dx) {

			return adjust_coefficients(f3_coefficients<Tnumber>, dx * dx);

		}

		// Forward difference; 4th order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> f4_coefficients{
			 15.0 / 4.0,
			-77.0 / 6.0,
			 107.0 / 6.0,
			-13.0,
			 61.0 / 12.0,
			-5.0 / 6.0
		};

		template<typename Tnumber>
		const std::vector<Tnumber> f4(const Tnumber dx) {

			return adjust_coefficients(f4_coefficients<Tnumber>, dx * dx);

		}

		// Backward difference; 1st order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> b1(const Tnumber dx) {

			return reverse_order(f1<Tnumber>(dx));

		}

		// Backward difference; 2nd order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> b2(const Tnumber dx) {

			return reverse_order(f2<Tnumber>(dx));

		}

		// Backward difference; 3rd order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> b3(const Tnumber dx) {

			return reverse_order(f3<Tnumber>(dx));

		}

		// Backward difference; 4th order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> b4(const Tnumber dx) {

			return reverse_order(f4<Tnumber>(dx));

		}

	}

	// Finite difference representation on non-uniform grid.
	// See Sundqvist and Veronis (1970).
	namespace nonuniform {

		// dx_vector; step size vector with four elements:
		// [0] dx_m2 = x_m1 - x_m2  TODO: Correct sign?
		// [1] dx_m1 = x - x_m1     TODO: Correct sign?
		// [2] dx_p1 = x_p1 - x  
		// [3] dx_p2 = x_p2 - x_p1

		// Central difference; ~2nd order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> c2(const std::vector<Tnumber>& dx_vector) {

			const Tnumber dx_m1 = dx_vector[1];
			const Tnumber dx_p1 = dx_vector[2];

			std::vector<Tnumber> row(3, 0.0);

			const Tnumber denominator = dx_p1 * dx_m1 * (1.0 + dx_p1 / dx_m1);

			// Sub-diagonal element.
			row[0] = 2.0 * (dx_p1 / dx_m1) / denominator;

			// Main diagonal element.
			row[1] = -2.0 * (1.0 + dx_p1 / dx_m1) / denominator;

			// Super-diagonal element.
			row[2] = 2.0 / denominator;

			return row;

		}

		// Central difference; ~4th order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> c4(const std::vector<Tnumber>& dx_vector) {

			const Tnumber dx_m2 = dx_vector[0];
			const Tnumber dx_m1 = dx_vector[1];
			const Tnumber dx_p1 = dx_vector[2];
			const Tnumber dx_p2 = dx_vector[3];

			std::vector<Tnumber> row(5, 0.0);

			const Tnumber denominator =
				-(dx_m1 + dx_m2) * pow(dx_p1 + dx_p2, 2) / 2.0
				+ 16.0 * pow(dx_p1, 2) * dx_m1
				+ 16.0 * dx_p1 * pow(dx_m1, 2)
				- pow(dx_m1 + dx_m2, 2) * (dx_p1 + dx_p2) / 2.0;

			// Coefficient of 2nd sub-diagonal.
			row[0] = -(dx_p1 + dx_p2) / denominator;

			// Coefficient of 1st sub-diagonal.
			row[1] = 32.0 * dx_p1 / denominator;

			// Coefficient of main diagonal.
			row[2] = -(-(dx_m1 + dx_m2) + 32.0 * dx_m1
				+ 32.0 * dx_p1 - (dx_p1 + dx_p2)) / denominator;

			// Coefficient of 1st super-diagonal.
			row[3] = 32.0 * dx_m1 / denominator;

			// Coefficient of 2nd super-diagonal.
			row[4] = -(dx_m1 + dx_m2) / denominator;

			return row;

		}

		// Forward difference; 1st order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> f1(const std::vector<Tnumber>& dx_vector) {

			const Tnumber dx_p1 = dx_vector[2];
			const Tnumber dx_p2 = dx_vector[3];

			std::vector<Tnumber> row(3, 0.0);

			const Tnumber denominator = dx_p1 * dx_p2 * (dx_p1 + dx_p2) / 2.0;

			// Coefficient of main diagonal.
			row[0] = dx_p2 / denominator;

			// Coefficient of 1st super-diagonal.
			row[1] = -(dx_p1 + dx_p2) / denominator;

			// Coefficient of 2nd super-diagonal.
			row[2] = dx_p1 / denominator;

			return row;

		}

		// Backward difference; 1st order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> b1(const std::vector<Tnumber>& dx_vector) {

			return reverse_order(f1<Tnumber>(dx_vector));

		}

	}

}

#if false
// Reverse order of coefficients and multiply by scalar.
template<typename Tnumber>
std::vector<Tnumber> reverse_orderTemplate(
	std::vector<Tnumber> coefficients,
	const Tnumber scalar) {

	std::reverse(coefficients.begin(), coefficients.end());

	for (auto& element : coefficients) {
		element *= scalar;
	}

	return coefficients;

}


// Divide coefficients by denominator. TODO: Change function name.
template<typename Tnumber>
std::vector<Tnumber> adjust_coefficientsTemplate(
	std::vector<Tnumber> coefficients,
	const Tnumber denominator) {

	for (auto& element : coefficients) {
		element /= denominator;
	}

	return coefficients;

}
#endif

// ###############################################################################


// Finite difference coefficients for first order derivative operator.
namespace coef_x1 {

	// Finite difference representation on uniform grid.
	// See Fornberg (1988).
	namespace uniform {

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

	// Finite difference representation on non-uniform grid.
	// See Sundqvist and Veronis (1970).
	namespace nonuniform {

		// dx_vector; step size vector with four elements:
		// [0] dx_m2 = x_m1 - x_m2  TODO: Correct sign?
		// [1] dx_m1 = x - x_m1     TODO: Correct sign?
		// [2] dx_p1 = x_p1 - x  
		// [3] dx_p2 = x_p2 - x_p1

		// Central difference; 2nd order accuracy.
		std::vector<double> c2(const std::vector<double>& dx_vector) {

			const double dx_m1 = dx_vector[1];
			const double dx_p1 = dx_vector[2];

			std::vector<double> row(3, 0.0);

			const double denominator = dx_p1 * (1.0 + dx_p1 / dx_m1);

			// Coefficient of 1st sub-diagonal.
			row[0] = -pow(dx_p1 / dx_m1, 2) / denominator;

			// Coefficient of main diagonal.
			row[1] = -(1.0 - pow(dx_p1 / dx_m1, 2)) / denominator;

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
			row[2] = 
				-(pow(dx_m1 + dx_m2, 2) - 32.0 * pow(dx_m1, 2) 
				  + 32.0 * pow(dx_p1, 2) - pow(dx_p1 + dx_p2, 2)) / denominator;

			// Coefficient of 1st super-diagonal.
			row[3] = -32.0 * pow(dx_m1, 2) / denominator;

			// Coefficient of 2nd super-diagonal.
			row[4] = pow(dx_m1 + dx_m2, 2) / denominator;

			return row;

		}

		// Forward difference; 1st order accuracy.
		std::vector<double> f1(const std::vector<double>& dx_vector) {

			const double dx_p1 = dx_vector[2];

			std::vector<double> row(2, 0.0);

			const double denominator = dx_p1;

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

			const double denominator = 
				pow(dx_p1, 2) * (dx_p1 + dx_p2) - dx_p1 * pow(dx_p1 + dx_p2, 2);

			// Coefficient of main diagonal.
			row[0] = (pow(dx_p2, 2) + 2 * dx_p1 * dx_p2) / denominator;

			// Coefficient of 1st super-diagonal.
			row[1] = -pow(dx_p1 + dx_p2, 2) / denominator;

			// Coefficient of 2nd super-diagonal.
			row[2] = pow(dx_p1, 2) / denominator;

			return row;

		}

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
namespace coef_x2 {

	// Finite difference representation on uniform grid.
	// See Fornberg (1988).
	namespace uniform {

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

	// Finite difference representation on non-uniform grid.
	// See Sundqvist and Veronis (1970).
	namespace nonuniform {

		// dx_vector; step size vector with four elements:
		// [0] dx_m2 = x_m1 - x_m2  TODO: Correct sign?
		// [1] dx_m1 = x - x_m1     TODO: Correct sign?
		// [2] dx_p1 = x_p1 - x  
		// [3] dx_p2 = x_p2 - x_p1

		// Central difference; ~2nd order accuracy.
		std::vector<double> c2(const std::vector<double>& dx_vector) {

			const double dx_m1 = dx_vector[1];
			const double dx_p1 = dx_vector[2];

			std::vector<double> row(3, 0.0);

			const double denominator = dx_p1 * dx_m1 * (1.0 + dx_p1 / dx_m1);

			// Sub-diagonal element.
			row[0] = 2.0 * (dx_p1 / dx_m1) / denominator;

			// Main diagonal element.
			row[1] = -2.0 * (1.0 + dx_p1 / dx_m1) / denominator;

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
			row[2] = -(-(dx_m1 + dx_m2) + 32.0 * dx_m1 
					   + 32.0 * dx_p1 - (dx_p1 + dx_p2)) / denominator;

			// Coefficient of 1st super-diagonal.
			row[3] = 32.0 * dx_m1 / denominator;

			// Coefficient of 2nd super-diagonal.
			row[4] = -(dx_m1 + dx_m2) / denominator;

			return row;

		}

		// Forward difference; 1st order accuracy.
		std::vector<double> f1(const std::vector<double>& dx_vector) {

			const double dx_p1 = dx_vector[2];
			const double dx_p2 = dx_vector[3];

			std::vector<double> row(3, 0.0);

			const double denominator = dx_p1 * dx_p2 * (dx_p1 + dx_p2) / 2.0;

			// Coefficient of main diagonal.
			row[0] = dx_p2 / denominator;

			// Coefficient of 1st super-diagonal.
			row[1] = -(dx_p1 + dx_p2) / denominator;

			// Coefficient of 2nd super-diagonal.
			row[2] = dx_p1 / denominator;

			return row;

		}

		// Backward difference; 1st order accuracy.
		std::vector<double> b1(const std::vector<double>& dx_vector) {

			return reverse_order(f1(dx_vector));

		}

	}

}


// Reverse order of coefficients and multiply by scalar.
std::vector<double> reverse_order(
	std::vector<double> coefficients, 
	const double scalar) {

	std::reverse(coefficients.begin(), coefficients.end());

	for (auto& element : coefficients) {
		element *= scalar;
	}

	return coefficients;

}


// Divide coefficients by denominator. TODO: Change function name.
std::vector<double> adjust_coefficients(
	std::vector<double> coefficients, 
	const double denominator) {

	for (auto& element : coefficients) {
		element /= denominator;
	}

	return coefficients;

}
