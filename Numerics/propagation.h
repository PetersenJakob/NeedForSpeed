#pragma once

#include <functional>
#include <stdexcept>
#include <vector>

#include "grid.h"
#include "propagator.h"

namespace propagation {

	namespace theta_1d {

		template <class T>
		void full(
			const std::vector<double>& time_grid,
			T& derivative,
			std::vector<double>& func,
			const double theta=0.5) {

			T identity = derivative.identity();

			for (int i = 0; i != time_grid.size() - 1; ++i) {

				double dt = time_grid[i + 1] - time_grid[i];

				propagator::theta_1d::full(dt, identity, derivative, func, theta);

			}

		}

	}

	namespace adi {

		template <class T1, class T2>
		void dr_2d(
			const std::vector<double>& time_grid,
			T1& derivative_1,
			T2& derivative_2,
			std::vector<double>& func,
			const double theta = 0.5) {

			T1 identity_1 = derivative_1.identity();
			T2 identity_2 = derivative_2.identity();

			for (int i = 0; i != time_grid.size() - 1; ++i) {

				double dt = time_grid[i + 1] - time_grid[i];

				propagator::adi::dr_2d(
					dt, 
					identity_1, identity_2,
					derivative_1, derivative_2,
					func, 
					theta);

			}

		}

		template <class T1, class T2>
		void cs_2d(
			const std::vector<double>& time_grid,
			T1& derivative_1,
			T2& derivative_2,
			MixedDerivative<T1, T2>& mixed,
			std::vector<double>& func,
			const double theta = 0.5,
			const double lambda = 0.5, 
			const int n_iterations = 1) {

			T1 identity_1 = derivative_1.identity();
			T2 identity_2 = derivative_2.identity();

			for (int i = 0; i != time_grid.size() - 1; ++i) {

				double dt = time_grid[i + 1] - time_grid[i];

				propagator::adi::cs_2d(
					dt,
					identity_1, identity_2,
					derivative_1, derivative_2,
					mixed
					func,
					theta, lambda, n_iterations);

			}

		}

		template <class T1, class T2, class T3>
		void dr_3d(
			const std::vector<double>& time_grid,
			T1& derivative_1,
			T2& derivative_2,
			T3& derivative_3,
			std::vector<double>& func,
			const double theta = 0.5) {

			T1 identity_1 = derivative_1.identity();
			T2 identity_2 = derivative_2.identity();
			T3 identity_3 = derivative_3.identity();

			for (int i = 0; i != time_grid.size() - 1; ++i) {

				double dt = time_grid[i + 1] - time_grid[i];

				propagator::adi::dr_3d(
					dt,
					identity_1, identity_2, identity_3,
					derivative_1, derivative_2, derivative_3,
					func,
					theta);

			}

		}

		template <class T>
		void dr_nd(
			const std::vector<double>& time_grid,
			std::vector<T>& derivative,
			std::vector<double>& func,
			const double theta = 0.5) {

			std::vector<T> identity = { derivative[0].identity() };

			for (int i = 1; i != identity.size(); ++i) {

				identity.push_back(derivative[i].identity());

			}

			for (int i = 0; i != time_grid.size() - 1; ++i) {

				double dt = time_grid[i + 1] - time_grid[i];

				propagator::adi::dr_nd(dt, identity, derivative, func, theta);

			}

		}

	}

}
