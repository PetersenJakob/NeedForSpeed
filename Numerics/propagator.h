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
		std::vector<double>& func,
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
		func = rhs * func;

		// Solve matrix equation.
		solver::band(lhs, func);

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
			std::vector<std::vector<double>>& func,
			const double theta_parameter = 0.5) {

			const int n_points_1 = identity_operator_1.order();
			const int n_points_2 = identity_operator_2.order();

			std::vector<double> inner(n_points_2, 0.0);
			std::vector<std::vector<double>> func_tmp(n_points_1, inner);

			// Eq. (2.68), right-hand-side.
			T2 lhs_2 = derivative_operator_2;
			lhs_2 *= time_step;
			for (int i = 0; i != n_points_1; ++i) {
				func_tmp[i] = lhs_2 * func[i];
			}

			T1 lhs_1 = derivative_operator_1;
			lhs_1 *= (1.0 - theta_parameter) * time_step;
			lhs_1 += identity_operator_1;

			for (int i = 0; i != n_points_2; ++i) {

				std::vector<double> func_strip(n_points_1, 0.0);

				for (int j = 0; j != n_points_1; ++i) {
					func_strip[j] = func[j][i];
				}

				func_strip = lhs_1 * func_strip;

				for (int j = 0; j != n_points_1; ++i) {
					func[j][i] = func_tmp[j][i] + func_strip[j];
				}

			}

			// Eq. (2.68), left-hand-side.
			T1 rhs_1 = derivative_operator_1;
			rhs_1 *= -theta_parameter * time_step;
			rhs_1 += identity_operator_1;

			for (int i = 0; i != n_points_2; ++i) {
			
				std::vector<double> func_strip(n_points_1, 0.0);

				for (int j = 0; j != n_points_1; ++i) {
					func_strip[j] = func[j][i];
				}

				solver::band(rhs_1, func_strip);

				for (int j = 0; j != n_points_1; ++i) {
					func[j][i] = func_strip[j];
				}

			}

			// Eq. (2.69), left-hand-side.
			for (int i = 0; i != n_points_1; ++i) {
				for (int j = 0; j != n_points_2; ++j) {
					func[i][j] -= theta_parameter * func_tmp[i][j];
				}
			}

			// Eq. (2.69), left-hand-side.
			T2 rhs_2 = derivative_operator_2;
			rhs_2 *= -theta_parameter;
			rhs_2 += identity_operator_2;

			for (int i = 0; i != n_points_1; ++i) {
				solver::band(rhs_2, func[i]);
			}

		}

		// Craig-Sneyd 2-dimensional scheme.
		template <class T1, class T2>
		void cs_2d(
			const double time_step,
			const T1& identity_operator_1,
			const T2& identity_operator_2,
			const T1& derivative_operator_1,
			const T2& derivative_operator_2,
			std::vector<std::vector<double>>& func,
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
			std::vector<std::vector<std::vector<double>>>& func,
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
			std::vector<std::vector<std::vector<double>>>& func,
			const double theta_parameter = 0.5) {

		}

		// Douglas-Rachford N-dimensional scheme.
		template <class T>
		void dr_nd(
			const double time_step,
			const std::vector<T>& identity_operator,
			const std::vector<T>& derivative_operator,
			std::vector<double>& func,
			const double theta_parameter = 0.5) {

			const int n_dimensions = (int)identity_operator.size();
			
			std::vector<int> n_points(n_dimensions, 0);
			for (int i = 0; i != n_dimensions; ++i) {
				n_points[i] = identity_operator[i].order();
			}

		}

		// Craig-Sneyd N-dimensional scheme.
		template <class T>
		void cs_nd(
			const double time_step,
			const std::vector<T>& identity_operator,
			const std::vector<T>& derivative_operator,
			std::vector<double>& func,
			const double theta_parameter = 0.5) {

			const int n_dimensions = (int)identity_operator.size();

			std::vector<int> n_points(n_dimensions, 0);
			for (int i = 0; i != n_dimensions; ++i) {
				n_points[i] = identity_operator[i].order();
			}

		}

	}

}
