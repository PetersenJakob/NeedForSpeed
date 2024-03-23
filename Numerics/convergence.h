#pragma once

#include <iostream>
#include <iomanip>

#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

#include "norm.h"
#include "propagation.h"


void grid_increment(
	const int n_increments,
	std::vector<double>& grid,
	std::function<std::vector<double>(double, double, int)> grid_generator);


std::vector<std::vector<double>> norm_vector(const int n_iterations);


double average_grid_spacing(const std::vector<double>& grid);


namespace convergence {

	template <class T>
	std::vector<std::vector<double>>
		theta_1d(
			std::vector<double>& time_grid,
			std::vector<std::vector<double>>& spatial_grid,

			std::function<std::vector<double>
				(const double, const double, const int)> grid_generator,

			std::function<T(std::vector<double>)> derivative_generator,

			std::function<std::vector<double>
				(double, std::vector<std::vector<double>>&)> solution_generator,

			std::string dimension,
			const int n_iterations,
			const int n_increments,
			const double theta = 0.5) {

		std::vector<std::vector<double>> norm = norm_vector(n_iterations);

		std::vector<double> func;
		std::vector<double> solution;

		for (int i = 0; i != n_iterations; ++i) {

			T derivative = derivative_generator(spatial_grid[0]);

			func = solution_generator(time_grid.front(), spatial_grid);

			propagation::theta_1d::full(time_grid, derivative, func, theta);

			solution = solution_generator(time_grid.back(), spatial_grid);

			std::vector<double> diff = norm::vector_diff(solution, func);

			if (dimension == "time") {
				norm[0][i] = average_grid_spacing(time_grid);
			}
			else if (dimension == "space") {
				norm[0][i] = average_grid_spacing(spatial_grid[0]);
			}
			else {
				throw std::invalid_argument("dimension unknown.");
			}

			norm[1][i] = norm::vector::infinity(diff);

			norm[2][i] = norm::vector::l1(diff);

			norm[3][i] = norm::vector::l2(diff);

			norm[4][i] = norm::function::l1(spatial_grid[0], diff);

			norm[5][i] = norm::function::l2(spatial_grid[0], diff);



			if (false) {
				std::cout << std::scientific << std::setprecision(5);
				std::cout << "Average time step: " << time_grid[1] - time_grid[0]
					<< "\t Time interval: " << (time_grid.back() - time_grid.front()) << std::endl;

				std::vector<double> initial_func = solution_generator(time_grid.front(), spatial_grid);

				for (int m = 0; m != spatial_grid[0].size(); ++m) {
					std::cout
						<< std::setw(3) << m
						<< std::setw(14) << spatial_grid[0][m]
						<< std::setw(14) << initial_func[m]
						<< std::setw(14) << solution[m]
						<< std::setw(14) << func[m]
						<< std::setw(14) << abs(solution[m] - func[m])
						<< std::endl;

				}
				std::cout << std::endl;
			}



			if (dimension == "time") {
				grid_increment(n_increments, time_grid, grid_generator);
			}
			else if (dimension == "space") {
				grid_increment(n_increments, spatial_grid[0], grid_generator);
			}
			else {
				throw std::invalid_argument("dimension unknown.");
			}

		}

		return norm;

	}

	namespace adi {

		template <class T1, class T2>
		std::vector<std::vector<double>>
			dr_2d(
				std::vector<double>& time_grid,
				std::vector<std::vector<double>>& spatial_grid,
				
				std::function<std::vector<double>
					(const double, const double, const int)> grid_generator,

				std::function<T1(std::vector<double>)> derivative_generator_1,
				std::function<T2(std::vector<double>)> derivative_generator_2,

				std::function<std::vector<double>
					(double, std::vector<std::vector<double>>&)> solution_generator,

				std::string dimension,
				const int n_iterations,
				const int n_increments,
				const double theta = 0.5) {

			std::vector<std::vector<double>> norm = norm_vector(n_iterations);

			std::vector<double> func;
			std::vector<double> solution;

			for (int i = 0; i != n_iterations; ++i) {

				T1 derivative_1 = derivative_generator_1(spatial_grid[0]);
				T2 derivative_2 = derivative_generator_2(spatial_grid[1]);

				func = solution_generator(time_grid.front(), spatial_grid);

				propagation::adi::dr_2d(time_grid,
					derivative_1, derivative_2,
					func, theta);

				solution = solution_generator(time_grid.back(), spatial_grid);

				std::vector<double> diff = norm::vector_diff(solution, func);

				if (dimension == "time") {
					norm[0][i] = average_grid_spacing(time_grid);
				}
				else if (dimension == "space_1") {
					norm[0][i] = average_grid_spacing(spatial_grid[0]);
				}
				else if (dimension == "space_2") {
					norm[0][i] = average_grid_spacing(spatial_grid[1]);
				}
				else {
					throw std::invalid_argument("dimension unknown.");
				}

				norm[1][i] = norm::vector::infinity(diff);

				norm[2][i] = norm::vector::l1(diff);

				norm[3][i] = norm::vector::l2(diff);

				norm[4][i] = norm::function::l1(spatial_grid[0], spatial_grid[1], diff);

				norm[5][i] = norm::function::l2(spatial_grid[0], spatial_grid[1], diff);

				if (dimension == "time") {
					grid_increment(n_increments, time_grid, grid_generator);
				}
				else if (dimension == "space_1") {
					grid_increment(n_increments, spatial_grid[0], grid_generator);
				}
				else if (dimension == "space_2") {
					grid_increment(n_increments, spatial_grid[1], grid_generator);
				}
				else {
					throw std::invalid_argument("dimension unknown.");
				}

			}

			return norm;

		}


		template <class T1, class T2>
		std::vector<std::vector<double>>
			dr_2d(
				std::vector<double>& time_grid,
				std::vector<std::vector<double>>& spatial_grid,

				std::function<std::vector<double>
				(const double, const double, const int)> grid_generator,

				const std::vector<double>& prefactors_1,
				const std::vector<double>& prefactors_2,

				std::vector< std::function<T1(std::vector<double>)> > derivative_generator_1,
				std::vector< std::function<T2(std::vector<double>)> > derivative_generator_2,

				std::function<std::vector<double>
				(double, std::vector<std::vector<double>>&)> solution_generator,

				std::string dimension,
				const int n_iterations,
				const int n_increments,
				const double theta = 0.5) {

			std::vector<std::vector<double>> norm = norm_vector(n_iterations);

			std::vector<double> func;
			std::vector<double> solution;

			for (int i = 0; i != n_iterations; ++i) {

				std::vector<T1> derivatives_1;
				derivatives_1.push_back(derivative_generator_1[0](spatial_grid[0]).identity());
				derivatives_1.push_back(derivative_generator_1[0](spatial_grid[0]));
				derivatives_1.push_back(derivative_generator_1[1](spatial_grid[0]));

				std::vector<T2> derivatives_2;
				derivatives_2.push_back(derivative_generator_2[0](spatial_grid[1]).identity());
				derivatives_2.push_back(derivative_generator_2[0](spatial_grid[1]));
				derivatives_2.push_back(derivative_generator_2[1](spatial_grid[1]));

				func = solution_generator(time_grid.front(), spatial_grid);

				propagation::adi::dr_2d(time_grid,
					prefactors_1, prefactors_2,
					derivatives_1, derivatives_2,
					func, theta);

				solution = solution_generator(time_grid.back(), spatial_grid);

				std::vector<double> diff = norm::vector_diff(solution, func);

				if (dimension == "time") {
					norm[0][i] = average_grid_spacing(time_grid);
				}
				else if (dimension == "space_1") {
					norm[0][i] = average_grid_spacing(spatial_grid[0]);
				}
				else if (dimension == "space_2") {
					norm[0][i] = average_grid_spacing(spatial_grid[1]);
				}
				else {
					throw std::invalid_argument("dimension unknown.");
				}

				norm[1][i] = norm::vector::infinity(diff);

				norm[2][i] = norm::vector::l1(diff);

				norm[3][i] = norm::vector::l2(diff);

				norm[4][i] = norm::function::l1(spatial_grid[0], spatial_grid[1], diff);

				norm[5][i] = norm::function::l2(spatial_grid[0], spatial_grid[1], diff);

				if (dimension == "time") {
					grid_increment(n_increments, time_grid, grid_generator);
				}
				else if (dimension == "space_1") {
					grid_increment(n_increments, spatial_grid[0], grid_generator);
				}
				else if (dimension == "space_2") {
					grid_increment(n_increments, spatial_grid[1], grid_generator);
				}
				else {
					throw std::invalid_argument("dimension unknown.");
				}

			}

			return norm;

		}


		template <class T1, class T2>
		std::vector<std::vector<double>>
			dr_2d(
				std::vector<double>& time_grid,
				std::vector<std::vector<double>>& spatial_grid,

				std::function<std::vector<double>
				(const double, const double, const int)> grid_generator,

				std::function<std::vector<std::vector<double>>
				(const std::vector<std::vector<double>>&)> prefactor_generator_1,
				std::function<std::vector<std::vector<double>>
				(const std::vector<std::vector<double>>&)> prefactor_generator_2,

				std::vector< std::function<T1(std::vector<double>)> > derivative_generator_1,
				std::vector< std::function<T2(std::vector<double>)> > derivative_generator_2,

				std::function<std::vector<double>
				(double, std::vector<std::vector<double>>&)> solution_generator,

				std::string dimension,
				const int n_iterations,
				const int n_increments,
				const double theta = 0.5) {

			std::vector<std::vector<double>> norm = norm_vector(n_iterations);

			std::vector<double> func;
			std::vector<double> solution;

			for (int i = 0; i != n_iterations; ++i) {

				std::vector<std::vector<double>> prefactors_1 = prefactor_generator_1(spatial_grid);
				std::vector<std::vector<double>> prefactors_2 = prefactor_generator_2(spatial_grid);

				std::vector<T1> derivatives_1;
				derivatives_1.push_back(derivative_generator_1[0](spatial_grid[0]).identity());
				derivatives_1.push_back(derivative_generator_1[0](spatial_grid[0]));
				derivatives_1.push_back(derivative_generator_1[1](spatial_grid[0]));

				std::vector<T2> derivatives_2;
				derivatives_2.push_back(derivative_generator_2[0](spatial_grid[1]).identity());
				derivatives_2.push_back(derivative_generator_2[0](spatial_grid[1]));
				derivatives_2.push_back(derivative_generator_2[1](spatial_grid[1]));

				func = solution_generator(time_grid.front(), spatial_grid);

				propagation::adi::dr_2d(time_grid,
					prefactors_1, prefactors_2,
					derivatives_1, derivatives_2,
					func, theta);

				solution = solution_generator(time_grid.back(), spatial_grid);

				std::vector<double> diff = norm::vector_diff(solution, func);

				if (dimension == "time") {
					norm[0][i] = average_grid_spacing(time_grid);
				}
				else if (dimension == "space_1") {
					norm[0][i] = average_grid_spacing(spatial_grid[0]);
				}
				else if (dimension == "space_2") {
					norm[0][i] = average_grid_spacing(spatial_grid[1]);
				}
				else {
					throw std::invalid_argument("dimension unknown.");
				}

				norm[1][i] = norm::vector::infinity(diff);

				norm[2][i] = norm::vector::l1(diff);

				norm[3][i] = norm::vector::l2(diff);

				norm[4][i] = norm::function::l1(spatial_grid[0], spatial_grid[1], diff);

				norm[5][i] = norm::function::l2(spatial_grid[0], spatial_grid[1], diff);

				if (dimension == "time") {
					grid_increment(n_increments, time_grid, grid_generator);
				}
				else if (dimension == "space_1") {
					grid_increment(n_increments, spatial_grid[0], grid_generator);
				}
				else if (dimension == "space_2") {
					grid_increment(n_increments, spatial_grid[1], grid_generator);
				}
				else {
					throw std::invalid_argument("dimension unknown.");
				}

			}

			return norm;

		}


		template <class T1, class T2>
		std::vector<std::vector<double>>
			cs_2d(
				std::vector<double>& time_grid,
				std::vector<std::vector<double>>& spatial_grid,

				std::function<std::vector<double>
				(const double, const double, const int)> grid_generator,

				const std::vector<double>& prefactors_1,
				const std::vector<double>& prefactors_2,

				std::vector< std::function<T1(std::vector<double>)> > derivative_generator_1,
				std::vector< std::function<T2(std::vector<double>)> > derivative_generator_2,

				std::function<std::vector<double>
				(double, std::vector<std::vector<double>>&)> solution_generator,

				std::string dimension,
				const int n_iterations,
				const int n_increments,
				const double theta = 0.5) {

			std::vector<std::vector<double>> norm = norm_vector(n_iterations);

			std::vector<double> func;
			std::vector<double> solution;

			for (int i = 0; i != n_iterations; ++i) {

				std::vector<T1> derivatives_1;
				derivatives_1.push_back(derivative_generator_1[0](spatial_grid[0]).identity());
				derivatives_1.push_back(derivative_generator_1[0](spatial_grid[0]));
				derivatives_1.push_back(derivative_generator_1[1](spatial_grid[0]));

				std::vector<T2> derivatives_2;
				derivatives_2.push_back(derivative_generator_2[0](spatial_grid[1]).identity());
				derivatives_2.push_back(derivative_generator_2[0](spatial_grid[1]));
				derivatives_2.push_back(derivative_generator_2[1](spatial_grid[1]));

				
				// TODO: Need prefactor-generator as cs_2d function parameter.
				MixedDerivative<T1, T2> mixed(
					derivative_generator_1[0](spatial_grid[0]),
					derivative_generator_2[0](spatial_grid[1]));
				mixed.set_prefactors(0.0);


				func = solution_generator(time_grid.front(), spatial_grid);

				propagation::adi::cs_2d(time_grid,
					derivatives_1[2], derivatives_2[2],
					mixed,
					func, theta);

				solution = solution_generator(time_grid.back(), spatial_grid);

				std::vector<double> diff = norm::vector_diff(solution, func);

				if (dimension == "time") {
					norm[0][i] = average_grid_spacing(time_grid);
				}
				else if (dimension == "space_1") {
					norm[0][i] = average_grid_spacing(spatial_grid[0]);
				}
				else if (dimension == "space_2") {
					norm[0][i] = average_grid_spacing(spatial_grid[1]);
				}
				else {
					throw std::invalid_argument("dimension unknown.");
				}

				norm[1][i] = norm::vector::infinity(diff);

				norm[2][i] = norm::vector::l1(diff);

				norm[3][i] = norm::vector::l2(diff);

				norm[4][i] = norm::function::l1(spatial_grid[0], spatial_grid[1], diff);

				norm[5][i] = norm::function::l2(spatial_grid[0], spatial_grid[1], diff);

				if (dimension == "time") {
					grid_increment(n_increments, time_grid, grid_generator);
				}
				else if (dimension == "space_1") {
					grid_increment(n_increments, spatial_grid[0], grid_generator);
				}
				else if (dimension == "space_2") {
					grid_increment(n_increments, spatial_grid[1], grid_generator);
				}
				else {
					throw std::invalid_argument("dimension unknown.");
				}

			}

			return norm;

		}

	}

}
