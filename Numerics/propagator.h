#pragma once

#include <functional>
#include <stdexcept>
#include <typeinfo>
#include <vector>

#include "band_diagonal_matrix.h"


// Time propagation schemes.
namespace propagator {

	// Theta scheme.
	// Andersen and Piterbarg, REF...
	template <class T>
	void theta_1d(
		const double dt,
		const T& identity,
		const T& derivative,
		std::vector<double>& func,
		const double theta = 0.5) {

		// Left-hand-side operator, Eq. ().
		T lhs = derivative;
		lhs *= - theta * dt;
		lhs += identity;

		// Right-hand-side operator, Eq. ().
		T rhs = derivative;
		rhs *= (1.0 - theta) * dt;
		rhs += identity;

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
			const double dt,
			const T1& identity_1,
			const T2& identity_2,
			const T1& derivative_1,
			const T2& derivative_2,
			std::vector<std::vector<double>>& func,
			const double theta = 0.5) {

			const int n_points_1 = identity_1.order();
			const int n_points_2 = identity_2.order();

			std::vector<double> inner(n_points_2, 0.0);
			std::vector<std::vector<double>> func_tmp(n_points_1, inner);

			// Eq. (2.68), right-hand-side.
			T2 lhs_2 = derivative_2;
			lhs_2 *= dt;
			for (int i = 0; i != n_points_1; ++i) {
				func_tmp[i] = lhs_2 * func[i];
			}

			T1 lhs_1 = derivative_1;
			lhs_1 *= (1.0 - theta) * dt;
			lhs_1 += identity_1;

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
			T1 rhs_1 = derivative_1;
			rhs_1 *= -theta * dt;
			rhs_1 += identity_1;

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
					func[i][j] -= theta * func_tmp[i][j];
				}
			}

			// Eq. (2.69), left-hand-side.
			T2 rhs_2 = derivative_2;
			rhs_2 *= -theta;
			rhs_2 += identity_2;

			for (int i = 0; i != n_points_1; ++i) {
				solver::band(rhs_2, func[i]);
			}

		}

		// Craig-Sneyd 2-dimensional scheme.
		template <class T1, class T2>
		void cs_2d(
			const double dt,
			const T1& identity_1,
			const T2& identity_2,
			const T1& derivative_1,
			const T2& derivative_2,

			std::function<std::vector<std::vector<double>>(
				T1&, 
				T2&,
				std::vector<std::vector<double>>)> derivative_12,

			const double theta = 0.5,
			const double lambda = 0.5) {

			const int n_points_1 = identity_1.order();
			const int n_points_2 = identity_2.order();

			std::vector<double> inner(n_points_2, 0.0);
			std::vector<std::vector<double>> func_tmp_1(n_points_1, inner);
			std::vector<std::vector<double>> func_tmp_2(n_points_1, inner);

			// ###############
			// Predictor step.
			// ###############

			// Eq. (2.68), right-hand-side.
			T2 lhs_2 = derivative_2;
			lhs_2 *= dt;
			for (int i = 0; i != n_points_1; ++i) {
				func_tmp_1[i] = lhs_2 * func[i];
			}


			// TODO: Should be d1dx1_1 and d1dx1_2!
			// TODO: What about prefix-function?
			func_tmp_2 = d2dxdy(derivative_1, derivative_2, func);


			T1 lhs_1 = derivative_1;
			lhs_1 *= (1.0 - theta) * dt;
			lhs_1 += identity_1;

			for (int i = 0; i != n_points_2; ++i) {

				std::vector<double> func_strip(n_points_1, 0.0);

				for (int j = 0; j != n_points_1; ++i) {
					func_strip[j] = func[j][i];
				}

				func_strip = lhs_1 * func_strip;

				for (int j = 0; j != n_points_1; ++i) {
					func[j][i] = func_tmp_1[j][i] + func_tmp_2[j][i] + func_strip[j];
				}

			}

			// Eq. (2.68), left-hand-side.
			T1 rhs_1 = derivative_1;
			rhs_1 *= -theta * dt;
			rhs_1 += identity_1;

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
					func[i][j] -= theta * func_tmp_1[i][j];
				}
			}

			// Eq. (2.69), left-hand-side.
			T2 rhs_2 = derivative_2;
			rhs_2 *= - theta;
			rhs_2 += identity_2;

			for (int i = 0; i != n_points_1; ++i) {
				solver::band(rhs_2, func[i]);
			}

			// ###############
			// Corrector step.
			// ###############

			 

		}

		// Douglas-Rachford 3-dimensional scheme.
		template <class T1, class T2, class T3>
		void dr_3d(
			const double dt,
			const T1& identity_1,
			const T2& identity_2,
			const T3& identity_3,
			const T1& derivative_1,
			const T2& derivative_2,
			const T3& derivative_3,
			std::vector<std::vector<std::vector<double>>>& func,
			const double theta = 0.5) {

		}

		// Craig-Sneyd 3-dimensional scheme.
		template <class T1, class T2, class T3>
		void cs_3d(
			const double dt,
			const T1& identity_1,
			const T2& identity_2,
			const T3& identity_3,
			const T1& derivative_1,
			const T2& derivative_2,
			const T3& derivative_3,
			std::vector<std::vector<std::vector<double>>>& func,
			const double theta = 0.5) {

		}

		// Douglas-Rachford N-dimensional scheme.
		template <class T>
		void dr_nd(
			const double dt,
			const std::vector<T>& identity,
			const std::vector<T>& derivative,
			std::vector<double>& func,
			const double theta = 0.5) {

			const int n_dimensions = (int)identity.size();
			
			std::vector<int> n_points(n_dimensions, 0);
			for (int i = 0; i != n_dimensions; ++i) {
				n_points[i] = identity[i].order();
			}

		}

		// Craig-Sneyd N-dimensional scheme.
		template <class T>
		void cs_nd(
			const double dt,
			const std::vector<T>& identity,
			const std::vector<T>& derivative,
			std::vector<double>& func,
			const double theta = 0.5) {

			const int n_dimensions = (int)identity.size();

			std::vector<int> n_points(n_dimensions, 0);
			for (int i = 0; i != n_dimensions; ++i) {
				n_points[i] = identity[i].order();
			}

		}

	}

}
