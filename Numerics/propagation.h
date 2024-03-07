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
					mixed,
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

#if false
// Move to .cpp file in UnitTest

void grid_increment(
	const int n_increments,
	std::vector<double>& grid,
	std::function<std::vector<double>(double, double, int)> grid_generator) {

	const double grid_min = grid.front();
	const double grid_max = grid.back();
	const int n_points = (int)grid.size() + n_increments;

	grid = grid_generator(grid_min, grid_max, n_points);

}


std::vector<std::vector<double>> norm_vector(const int n_iterations) {

	std::vector<double> inner(n_iterations, 0.0);

	//
	// step_size_vec
	// max_norm
	// l1_vec_norm
	// l2_vec_norm
	// l1_func_norm
	// l2_func_norm
	//

	std::vector<std::vector<double>> outer(6, inner);

	return outer;

}


double average_grid_spacing(std::vector<double>& grid) {

	return (grid.back() - grid.front()) / (grid.size() - 1);

}


namespace convergence {

	template <class T>
	std::vector<std::vector<double>>
		theta_1d(
			std::vector<double>& time_grid,
			std::vector<double>& spatial_grid, 
			std::function<std::vector<double>(double, double, int)> grid_generator,
			
			// Make changes to all operators in derivatves.cpp, such that only argument is spatial grid
			std::function<T(std::vector<double>)> derivative_generator,
			
			std::vector<double>& initial_func,
			std::vector<double>& solution,
			std::string dimension,
			const int n_iterations,
			const int n_increments,
			const double theta = 0.5) {
	
		std::vector<std::vector<double>> norm = norm_vector(n_iterations);

		for (int i = 0; i != n_iterations; ++i) {

			if (dimension == "time") {
				grid_increment(n_increments, time_grid, grid_generator);
				norm[0][i] = average_grid_spacing(time_grid);
			}
			else if (dimension == "space") {
				grid_increment(n_increments, spatial_grid, grid_generator);
				norm[0][i] = average_grid_spacing(spatial_grid);
			}
			else {
				throw std::invalid_argument("dimension unknown.");
			}

			T derivative = derivative_generator(spatial_grid);

			std::vector<double> func = initial_func;
			propagation::theta_1d::full(time_grid, derivative, func, theta);

			std::vector<double> diff = norm::vector_diff(solution, func);

			norm[1][i] = norm::vector::infinity(diff);

			norm[2][i] = norm::vector::l1(diff);

			norm[3][i] = norm::vector::l2(diff);

			norm[4][i] = norm::function::l1(spatial_grid, diff);

			norm[5][i] = norm::function::l2(spatial_grid, diff);

		}

		return norm;

	}

	namespace adi {

		template <class T1, class T2>
		std::vector<std::vector<double>>
			dr_2d(
				std::vector<double>& time_grid,
				std::vector<std::vector<double>>& spatial_grid,
				std::function<std::vector<double>(double, double, int)> grid_generator,

				// Make changes to all operators in derivatves.cpp, such that only argument is spatial grid
				std::function<T1(std::vector<double>)> derivative_generator_1,
				std::function<T2(std::vector<double>)> derivative_generator_2,

				std::vector<double>& initial_func,
				std::vector<double>& solution,
				std::string dimension,
				const int n_iterations,
				const int n_increments) {

			std::vector<std::vector<double>> norm = norm_vector(n_iterations);

			for (int i = 0; i != n_iterations; ++i) {

				if (dimension == "time") {
					grid_increment(n_increments, time_grid, grid_generator);
					norm[0][i] = average_grid_spacing(time_grid);
				}
				else if (dimension == "space_1") {
					grid_increment(n_increments, spatial_grid[0], grid_generator);
					norm[0][i] = average_grid_spacing(spatial_grid[0]);
				}
				else if (dimension == "space_2") {
					grid_increment(n_increments, spatial_grid[1], grid_generator);
					norm[0][i] = average_grid_spacing(spatial_grid[1]);
				}
				else {
					throw std::invalid_argument("dimension unknown.");
				}

				T1 derivative_1 = derivative_generator_1(spatial_grid[0]);
				T2 derivative_2 = derivative_generator_2(spatial_grid[1]);

				std::vector<double> func = initial_func;
				propagation::adi::dr_2d(time_grid, 
					derivative_1, derivative_2, 
					func, theta);

				std::vector<double> diff = norm::vector_diff(solution, func);

				norm[1][i] = norm::vector::infinity(diff);

				norm[2][i] = norm::vector::l1(diff);

				norm[3][i] = norm::vector::l2(diff);

				norm[4][i] = norm::function::l1(spatial_grid[0], diff);

				norm[5][i] = norm::function::l2(spatial_grid[0], diff);

			}

			return norm;

		}

	}

}
#endif