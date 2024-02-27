#pragma once

#include <functional>
#include <stdexcept>
#include <typeinfo>
#include <vector>

#include "band_diagonal_matrix.h"


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

			// Solve matrix equation.
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

		// Douglas-Rachford 2-dimensional scheme.
		// References
		// - AP: Andersen and Piterbarg (2010).
		template <class T1, class T2>
		void dr_2d(
			const double dt,
			const T1& identity_1,
			const T2& identity_2,
			const T1& derivative_1,
			const T2& derivative_2,
			std::vector<std::vector<double>>& func,
			const double theta = 0.5) {

			// TODO: Change func to ordinary std::vector<double>?

			const int n_points_1 = identity_1.order();
			const int n_points_2 = identity_2.order();

			std::vector<double> inner(n_points_2, 0.0);
			std::vector<std::vector<double>> func_tmp(n_points_1, inner);

			std::vector<double> func_strip(n_points_1, 0.0);

			// AP Eq. (2.68), right-hand-side.
			T2 rhs_2 = derivative_2;
			rhs_2 *= dt;

			for (int i = 0; i != n_points_1; ++i) {
				func_tmp[i] = rhs_2 * func[i];
			}

			T1 rhs_1 = derivative_1;
			rhs_1 *= (1.0 - theta) * dt;
			rhs_1 += identity_1;

			for (int i = 0; i != n_points_2; ++i) {

				std::vector<double> func_strip(n_points_1, 0.0);

				for (int j = 0; j != n_points_1; ++i) {
					func_strip[j] = func[j][i];
				}

				func_strip = rhs_1 * func_strip;

				for (int j = 0; j != n_points_1; ++i) {
					func[j][i] = func_strip[j] + func_tmp[j][i];
				}

			}

			// AP Eq. (2.68), left-hand-side.
			T1 lhs_1 = derivative_1;
			lhs_1 *= - theta * dt;
			lhs_1 += identity_1;

			for (int i = 0; i != n_points_2; ++i) {
			
				for (int j = 0; j != n_points_1; ++i) {
					func_strip[j] = func[j][i];
				}

				solver::band(lhs_1, func_strip);

				for (int j = 0; j != n_points_1; ++i) {
					func[j][i] = func_strip[j];
				}

			}

			// AP Eq. (2.69), right-hand-side.
			for (int i = 0; i != n_points_1; ++i) {
				for (int j = 0; j != n_points_2; ++j) {
					func[i][j] -= theta * func_tmp[i][j];
				}
			}

			// AP Eq. (2.69), left-hand-side.
			T2 lhs_2 = derivative_2;
			lhs_2 *= - theta * dt;
			lhs_2 += identity_2;

			for (int i = 0; i != n_points_1; ++i) {
				solver::band(lhs_2, func[i]);
			}

		}

		// Craig-Sneyd 2-dimensional scheme.
		// References
		// - AP: Andersen and Piterbarg (2010).
		template <class T1, class T2>
		void cs_2d(
			const double dt,
			const T1& identity_1,
			const T2& identity_2,
			const T1& derivative_1,
			const T2& derivative_2,

			std::function<std::vector<std::vector<double>>(
				std::vector<std::vector<double>>)> derivative_12,

			std::vector<std::vector<double>>& func,

			const double theta = 0.5,
			const double lambda = 0.5,
			const int n_iterations = 1) {

			const int n_points_1 = identity_1.order();
			const int n_points_2 = identity_2.order();

			std::vector<double> inner(n_points_2, 0.0);
			std::vector<std::vector<double>> func_tmp_1(n_points_1, inner);
			std::vector<std::vector<double>> func_tmp_2(n_points_1, inner);
			std::vector<std::vector<double>> func_tmp_3(n_points_1, inner);

			std::vector<double> func_strip(n_points_1, 0.0);

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

			for (int n = 0; n != n_interations; ++n) {

				// ###############
				// Predictor step.
				// ###############

				// AP Eq. (2.88), right-hand-side.
				for (int i = 0; i != n_points_2; ++i) {

					for (int j = 0; j != n_points_1; ++j) {
						func_strip[j] = func[j][i];
					}

					func_strip = rhs_1 * func_strip;

					for (int j = 0; j != n_points_1; ++j) {
						func_tmp_1[j][i] = func_strip[j];
					}

				}

				for (int i = 0; i != n_points_1; ++i) {
					func_tmp_2[i] = rhs_2 * func[i];
				}

				func_tmp_3 = derivative_12(func);

				for (int i = 0; i != n_points_1; ++i) {
					for (int j = 0; j != n_points_2; ++j) {
						func[i][j] = func_tmp_1[i][j] + func_tmp_2[i][j] 
							+ dt * func_tmp_3[i][j];
					}
				}

				// AP Eq. (2.88), left-hand-side.
				for (int i = 0; i != n_points_2; ++i) {

					for (int j = 0; j != n_points_1; ++j) {
						func_strip[j] = func[j][i];
					}

					solver::band(lhs_1, func_strip);

					for (int j = 0; j != n_points_1; ++j) {
						func[j][i] = func_strip[j];
					}

				}

				// AP Eq. (2.89), right-hand-side.
				for (int i = 0; i != n_points_1; ++i) {
					for (int j = 0; j != n_points_2; ++j) {
						func[i][j] -= theta * func_tmp_2[i][j];
					}
				}

				// AP Eq. (2.89), left-hand-side.
				for (int i = 0; i != n_points_1; ++i) {
					solver::band(lhs_2, func[i]);
				}

				// ###############
				// Corrector step.
				// ###############

				// AP Eq. (2.90), right-hand-side.
				func = derivative_12(func);

				for (int i = 0; i != n_points_1; ++i) {
					for (int j = 0; j != n_points_2; ++j) {
						func[i][j] *= lambda * dt;
						func[i][j] = func_tmp_1[i][j] + func_tmp_2[i][j]
							+ (1.0 - lambda) * dt * func_tmp_3[i][j];
					}
				}

				// AP Eq. (2.90), left-hand-side.
				for (int i = 0; i != n_points_2; ++i) {

					for (int j = 0; j != n_points_1; ++j) {
						func_strip[j] = func[j][i];
					}

					solver::band(lhs_1, func_strip);

					for (int j = 0; j != n_points_1; ++j) {
						func[j][i] = func_strip[j];
					}

				}

				// AP Eq. (2.91), right-hand-side.
				for (int i = 0; i != n_points_1; ++i) {
					for (int j = 0; j != n_points_2; ++j) {
						func[i][j] -= theta * func_tmp_2[i][j];
					}
				}

				// AP Eq. (2.91), left-hand-side.
				for (int i = 0; i != n_points_1; ++i) {
					solver::band(lhs_2, func[i]);
				}

			}

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
