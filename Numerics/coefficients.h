#pragma once

#include <vector>


// Finite difference coefficients for first order derivative operator.
namespace coef_x1Template {

	// Finite difference representation on uniform grid.
	namespace uniform {

		// Central difference; 2nd order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> c2(const Tnumber dx);

		// Central difference; 4th order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> c4(const Tnumber dx);

		// Forward difference; 1st order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> f1(const Tnumber dx);

		// Forward difference; 2nd order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> f2(const Tnumber dx);

		// Forward difference; 3rd order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> f3(const Tnumber dx);

		// Forward difference; 4th order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> f4(const Tnumber dx);

		// Backward difference; 1st order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> b1(const Tnumber dx);

		// Backward difference; 2nd order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> b2(const Tnumber dx);

		// Backward difference; 3rd order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> b3(const Tnumber dx);

		// Backward difference; 4th order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> b4(const Tnumber dx);

	}

	// Finite difference representation on non-uniform grid.
	namespace nonuniform {

		// Central difference; 2nd order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> c2(const std::vector<Tnumber>& dx_vector);

		// Central difference; 4th order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> c4(const std::vector<Tnumber>& dx_vector);

		// Forward difference; 1st order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> f1(const std::vector<Tnumber>& dx_vector);

		// Forward difference; 2nd order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> f2(const std::vector<Tnumber>& dx_vector);

		// Backward difference; 1st order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> b1(const std::vector<Tnumber>& dx_vector);

		// Backward difference; 2nd order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> b2(const std::vector<Tnumber>& dx_vector);

	}

}


// Finite difference coefficients for second order derivative operator.
namespace coef_x2Template {

	// Finite difference representation on uniform grid.
	namespace uniform {

		// Central difference; 2nd order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> c2(const Tnumber dx);

		// Central difference; 4th order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> c4(const Tnumber dx);

		// Forward difference; 1st order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> f1(const Tnumber dx);

		// Forward difference; 2nd order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> f2(const Tnumber dx);

		// Forward difference; 3rd order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> f3(const Tnumber dx);

		// Forward difference; 4th order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> f4(const Tnumber dx);

		// Backward difference; 1st order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> b1(const Tnumber dx);

		// Backward difference; 2nd order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> b2(const Tnumber dx);

		// Backward difference; 3rd order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> b3(const Tnumber dx);

		// Backward difference; 4th order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> b4(const Tnumber dx);

	}

	// Finite difference representation on non-uniform grid.
	namespace nonuniform {

		// Central difference; ~2nd order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> c2(const std::vector<Tnumber>& dx_vector);

		// Central difference; ~4th order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> c4(const std::vector<Tnumber>& dx_vector);

		// Forward difference; 1st order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> f1(const std::vector<Tnumber>& dx_vector);

		// TODO: Forward difference; ~2nd order accuracy.

		// Backward difference; 1st order accuracy.
		template<typename Tnumber>
		const std::vector<Tnumber> b1(const std::vector<Tnumber>& dx_vector);

		// TODO: Backward difference; ~2nd order accuracy.

	}

}


// Reverse order of coefficients and multiply by scalar.
template<typename Tnumber>
std::vector<Tnumber> reverse_orderTemplate(
	std::vector<Tnumber> coef,
	const Tnumber scalar = 1.0);


// Divide coefficients by denominator.
template<typename Tnumber>
std::vector<Tnumber> adjust_coefficientsTemplate(
	std::vector<Tnumber> coefficients,
	const Tnumber denominator);


// ###############################################################################


// Finite difference coefficients for first order derivative operator.
namespace coef_x1 {

	// Finite difference representation on uniform grid.
	namespace uniform {

		// Central difference; 2nd order accuracy.
		std::vector<double> c2(const double dx);

		// Central difference; 4th order accuracy.
		std::vector<double> c4(const double dx);

		// Forward difference; 1st order accuracy.
		std::vector<double> f1(const double dx);

		// Forward difference; 2nd order accuracy.
		std::vector<double> f2(const double dx);

		// Forward difference; 3rd order accuracy.
		std::vector<double> f3(const double dx);

		// Forward difference; 4th order accuracy.
		std::vector<double> f4(const double dx);

		// Backward difference; 1st order accuracy.
		std::vector<double> b1(const double dx);

		// Backward difference; 2nd order accuracy.
		std::vector<double> b2(const double dx);

		// Backward difference; 3rd order accuracy.
		std::vector<double> b3(const double dx);

		// Backward difference; 4th order accuracy.
		std::vector<double> b4(const double dx);

	}

	// Finite difference representation on non-uniform grid.
	namespace nonuniform {

		// Central difference; 2nd order accuracy.
		std::vector<double> c2(const std::vector<double>& dx_vector);

		// Central difference; 4th order accuracy.
		std::vector<double> c4(const std::vector<double>& dx_vector);

		// Forward difference; 1st order accuracy.
		std::vector<double> f1(const std::vector<double>& dx_vector);

		// Forward difference; 2nd order accuracy.
		std::vector<double> f2(const std::vector<double>& dx_vector);

		// Backward difference; 1st order accuracy.
		std::vector<double> b1(const std::vector<double>& dx_vector);

		// Backward difference; 2nd order accuracy.
		std::vector<double> b2(const std::vector<double>& dx_vector);

	}

}


// Finite difference coefficients for second order derivative operator.
namespace coef_x2 {

	// Finite difference representation on uniform grid.
	namespace uniform {

		// Central difference; 2nd order accuracy.
		std::vector<double> c2(const double dx);

		// Central difference; 4th order accuracy.
		std::vector<double> c4(const double dx);

		// Forward difference; 1st order accuracy.
		std::vector<double> f1(const double dx);

		// Forward difference; 2nd order accuracy.
		std::vector<double> f2(const double dx);

		// Forward difference; 3rd order accuracy.
		std::vector<double> f3(const double dx);

		// Forward difference; 4th order accuracy.
		std::vector<double> f4(const double dx);

		// Backward difference; 1st order accuracy.
		std::vector<double> b1(const double dx);

		// Backward difference; 2nd order accuracy.
		std::vector<double> b2(const double dx);

		// Backward difference; 3rd order accuracy.
		std::vector<double> b3(const double dx);

		// Backward difference; 4th order accuracy.
		std::vector<double> b4(const double dx);

	}

	// Finite difference representation on non-uniform grid.
	namespace nonuniform {

		// Central difference; ~2nd order accuracy.
		std::vector<double> c2(const std::vector<double>& dx_vector);

		// Central difference; ~4th order accuracy.
		std::vector<double> c4(const std::vector<double>& dx_vector);

		// Forward difference; 1st order accuracy.
		std::vector<double> f1(const std::vector<double>& dx_vector);

		// TODO: Forward difference; ~2nd order accuracy.

		// Backward difference; 1st order accuracy.
		std::vector<double> b1(const std::vector<double>& dx_vector);

		// TODO: Backward difference; ~2nd order accuracy.

	}

}


// Reverse order of coefficients and multiply by scalar.
std::vector<double> reverse_order(
	std::vector<double> coef, 
	const double scalar = 1.0);


// Divide coefficients by denominator.
std::vector<double> adjust_coefficients(
	std::vector<double> coefficients, 
	const double denominator);
