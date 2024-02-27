#pragma once

#include <functional>
#include <stdexcept>
#include <typeinfo>
#include <vector>

#include "band_diagonal_matrix.h"
#include "matrix_equation_solver.h"
#include "utility.h"


// Time propagation schemes.
namespace propagator {

	// Theta scheme.
	// References
	// - AP: Andersen and Piterbarg (2010).
	namespace theta_1d {

		template <class T>
		void step_1(
			const double dt,
			const T& identity,
			const T& derivative,
			std::vector<double>& func,
			const double theta = 0.5) {

			// AP Eq. (2.18), right-hand-side operator.
			// Operator evaluated at time t + dt (see AP Remarks 2.2.4).
			T rhs = derivative;
			rhs *= (1.0 - theta) * dt;
			rhs += identity;

			func = rhs * func;

		}

		template <class T>
		void step_2(
			const double dt,
			const T& identity,
			const T& derivative,
			std::vector<double>& func,
			const double theta = 0.5) {

			// AP Eq. (2.18), right-hand-side operator.
			// Operator evaluated at time t (see AP Remarks 2.2.4).
			T lhs = derivative;
			lhs *= -theta * dt;
			lhs += identity;

			solver::band(lhs, func);

		}

		template <class T>
		void full(
			const double dt,
			const T& identity,
			const T& derivative,
			std::vector<double>& func,
			const double theta = 0.5) {

			// Step one is carried out at time t + dt.
			step_1(dt, identity, derivative, func, theta);

			// Step one is carried out at time t.
			step_2(dt, identity, derivative, func, theta);

		}

	}

	// Alternating direction implicit schemes.
	namespace adi {

		// Douglas-Rachford scheme, 2-dimensional.
		// References
		// - AP: Andersen and Piterbarg (2010).
		template <class T1, class T2>
		void dr_2d(
			const double dt,
			const T1& identity_1,
			const T2& identity_2,
			const T1& derivative_1,
			const T2& derivative_2,
			std::vector<double>& func,
			const double theta = 0.5) {

			const int n_p_1 = identity_1.order();
			const int n_p_2 = identity_2.order();
			const int n_points = n_p_1 * n_p_2;

			std::vector<double> func_tmp_1(n_points, 0.0);
			std::vector<double> func_tmp_2(n_points, 0.0);

			// ##########
			// Operators.
			// ##########

			// AP Eq. (2.68) and (2.69), right-hand-side.
			T1 rhs_1 = derivative_1;
			rhs_1 *= (1.0 - theta) * dt;
			rhs_1 += identity_1;

			T2 rhs_2 = derivative_2;
			rhs_2 *= dt;

			// AP Eq. (2.68) and (2.69), left-hand-side.
			T1 lhs_1 = derivative_1;
			lhs_1 *= -theta * dt;
			lhs_1 += identity_1;

			T2 lhs_2 = derivative_2;
			lhs_2 *= -theta * dt;
			lhs_2 += identity_2;

			// ############
			// Propagation.
			// ############
			
			// AP Eq. (2.68), right-hand-side.
			func_tmp_1 = action_2d(n_p_1, n_p_2, 1, false, rhs_1, func);
			func_tmp_2 = action_2d(n_p_2, n_p_1, 2, false, rhs_2, func);

			for (int i = 0; i != n_points; ++i) {
				func[i] = func_tmp_1[i] + func_tmp_2[i];
			}

			// AP Eq. (2.68), left-hand-side.
			func = action_2d(n_p_1, n_p_2, 1, true, lhs_1, func);

			// AP Eq. (2.69), right-hand-side.
			for (int i = 0; i != n_points; ++i) {
				func[i] -= theta * func_tmp_2[i];
			}

			// AP Eq. (2.69), left-hand-side.
			func = action_2d(n_p_2, n_p_1, 2, true, lhs_2, func);

		}

		// Craig-Sneyd scheme, 2-dimensional.
		// References
		// - AP: Andersen and Piterbarg (2010).
		template <class T1, class T2>
		void cs_2d(
			const double dt,
			const T1& identity_1,
			const T2& identity_2,
			const T1& derivative_1,
			const T2& derivative_2,
			MixedDerivative<T1, T2>& mixed,
			std::vector<double>& func,
			const double theta = 0.5,
			const double lambda = 0.5,
			const int n_iterations = 1) {

			const int n_p_1 = identity_1.order();
			const int n_p_2 = identity_2.order();

			const int n_points = n_p_1 * n_p_2;

			std::vector<double> func_tmp_1(n_points, 0.0);
			std::vector<double> func_tmp_2(n_points, 0.0);
			std::vector<double> func_tmp_3(n_points, 0.0);

			// ##########
			// Operators.
			// ##########

			// AP Eq. (2.88) and (2.90), right-hand-side.
			T1 rhs_1 = derivative_1;
			rhs_1 *= (1.0 - theta) * dt;
			rhs_1 += identity_1;

			T2 rhs_2 = derivative_2;
			rhs_1 *= dt;

			// AP Eq. (2.88) and (2.90), left-hand-side.
			T1 lhs_1 = derivative_1;
			lhs_1 *= - theta * dt;
			lhs_1 += identity_1;

			T2 lhs_2 = derivative_2;
			lhs_2 *= - theta * dt;
			lhs_2 += identity_2;

			// ############
			// Propagation.
			// ############

			for (int n = 0; n != n_interations; ++n) {

				// ###############
				// Predictor step.
				// ###############

				// AP Eq. (2.88), right-hand-side.
				func_tmp_1 = action_2d(n_p_1, n_p_2, 1, false, rhs_1, func);
				func_tmp_2 = action_2d(n_p_2, n_p_1, 2, false, rhs_2, func);
				func_tmp_3 = mixed.d2dxdy(func);

				for (int i = 0; i != n_points; ++i) {
					func[i] = func_tmp_1[i] + func_tmp_2[i] + dt * func_tmp_3[i];
				}

				// AP Eq. (2.88), left-hand-side.
				func = action_2d(n_p_1, n_p_2, 1, true, lhs_1, func);

				// AP Eq. (2.89), right-hand-side.
				for (int i = 0; i != n_points; ++i) {
					func[i] -= theta * func_tmp_2[i];
				}

				// AP Eq. (2.89), left-hand-side.
				func = action_2d(n_p_2, n_p_1, 2, true, lhs_2, func);

				// ###############
				// Corrector step.
				// ###############

				// AP Eq. (2.90), right-hand-side.
				func = mixed(func);

				for (int i = 0; i != n_points; ++i) {
					func[i] *= lambda * dt;
					func[i] += func_tmp_1[i] + func_tmp_2[i] 
						+ (1.0 - lambda) * dt * func_tmp_3[i];
				}

				// AP Eq. (2.90), left-hand-side.
				func = action_2d(n_p_1, n_p_2, 1, true, lhs_1, func);

				// AP Eq. (2.91), right-hand-side.
				for (int i = 0; i != n_points; ++i) {
					func[i] -= theta * func_tmp_2[i];
				}

				// AP Eq. (2.91), left-hand-side.
				func = action_2d(n_p_2, n_p_1, 2, true, lhs_2, func);

			}

		}

		// Douglas-Rachford scheme, 3-dimensional.
		template <class T1, class T2, class T3>
		void dr_3d(
			const double dt,
			const T1& identity_1,
			const T2& identity_2,
			const T3& identity_3,
			const T1& derivative_1,
			const T2& derivative_2,
			const T3& derivative_3,
			std::vector<double>& func,
			const double theta = 0.5) {

		}

		// Craig-Sneyd scheme, 3-dimensional.
		template <class T1, class T2, class T3>
		void cs_3d(
			const double dt,
			const T1& identity_1,
			const T2& identity_2,
			const T3& identity_3,
			const T1& derivative_1,
			const T2& derivative_2,
			const T3& derivative_3,
			MixedDerivative<T1, T2>& mixed_12,
			MixedDerivative<T1, T3>& mixed_13,
			MixedDerivative<T2, T3>& mixed_23,
			std::vector<double>& func,
			const double theta = 0.5,
			const double lambda = 0.5,
			const int n_iterations = 1) {

		}

		// Douglas-Rachford scheme, N-dimensional.
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

		// Craig-Sneyd scheme, N-dimensional.
		template <class T>
		void cs_nd(
			const double dt,
			const std::vector<T>& identity,
			const std::vector<T>& derivative,
			std::vector<MixedDerivative<T, T>>& mixed,
			std::vector<double>& func,
			const double theta = 0.5,
			const double lambda = 0.5,
			const int n_iterations = 1) {

			const int n_dimensions = (int)identity.size();

			std::vector<int> n_points(n_dimensions, 0);
			for (int i = 0; i != n_dimensions; ++i) {
				n_points[i] = identity[i].order();
			}

		}

	}

}
