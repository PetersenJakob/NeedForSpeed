#pragma once

#include <stdexcept>
#include <typeinfo>
#include <vector>

#include "band_diagonal_matrix.h"


// Time propagation schemes.
namespace propagator {

	// Theta scheme.
	template <class T>
	void theta(
		const double time_step,
		const T& identity_operator,
		const T& derivative_operator,
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

	// Alternating direction implicit schemes.
	namespace adi {

		// Douglas-Rachford 2-dimensional scheme.
		template <class T1, class T2>
		void dr_2d(
			const double time_step,
			const T1& identity_operator_1,
			const T2& identity_operator_2,
			const T1& derivative_operator_1,
			const T2& derivative_operator_2,
			std::vector<std::vector<double>>& column,
			const double theta_parameter = 0.5) {

		}

		// Craig-Sneyd 2-dimensional scheme.
		template <class T1, class T2>
		void cs_2d(
			const double time_step,
			const T1& identity_operator_1,
			const T2& identity_operator_2,
			const T1& derivative_operator_1,
			const T2& derivative_operator_2,
			std::vector<std::vector<double>>& column,
			const double theta_parameter = 0.5) {

		}

		// Douglas-Rachford 3-dimensional scheme.
		template <class T1, class T2, class T3>
		void dr_3d(
			const double time_step,
			const T1& identity_operator_1,
			const T2& identity_operator_2,
			const T3& identity_operator_3,
			const T1& derivative_operator_1,
			const T2& derivative_operator_2,
			const T3& derivative_operator_3,
			std::vector<std::vector<std::vector<double>>>& column,
			const double theta_parameter = 0.5) {

		}

		// Craig-Sneyd 3-dimensional scheme.
		template <class T1, class T2, class T3>
		void cs_3d(
			const double time_step,
			const T1& identity_operator_1,
			const T2& identity_operator_2,
			const T3& identity_operator_3,
			const T1& derivative_operator_1,
			const T2& derivative_operator_2,
			const T3& derivative_operator_3,
			std::vector<std::vector<std::vector<double>>>& column,
			const double theta_parameter = 0.5) {

		}

		// Douglas-Rachford N-dimensional scheme.
		template <class T>
		void dr_2d(
			const double time_step,
			const std::vector<T>& identity_operator,
			const std::vector<T>& derivative_operator,
			std::vector<double>& column,
			const double theta_parameter = 0.5) {

		}

		// Craig-Sneyd N-dimensional scheme.
		template <class T>
		void cs_2d(
			const double time_step,
			const std::vector<T>& identity_operator,
			const std::vector<T>& derivative_operator,
			std::vector<double>& column,
			const double theta_parameter = 0.5) {

		}

	}

}
