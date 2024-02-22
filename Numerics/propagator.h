#pragma once


#include <stdexcept>
#include <typeinfo>
#include <vector>

#include "band_diagonal_matrix.h"

namespace propagator {


	template <class T>
	void theta(
		const double time_step,
		T& identity_operator,
		T& derivative_operator,
		std::vector<double>& column,
		const double theta_parameter = 0.5) {

		// Left-hand-side operator.
		T lhs = derivative_operator;
		lhs *= -theta_parameter * time_step;
		lhs += identity_operator;

		// Right-hand-side operator.
		T rhs = derivative_operator;
		rhs *= (1.0 - theta_parameter) * time_step;
		rhs += identity_operator;

		// Evaluation of righ-hand-side.
		column = rhs * column;

		// Solve matrix equation.
		solver::band(lhs, column);

	}

}