#include "pch.h"


std::vector<double> linear_regression(
	std::vector<std::vector<double>>& norm,
	const bool show_output = true,
	const std::string file_name = "slr.csv") {

	std::vector<double> inner(norm[0].size(), 0.0);
	std::vector<std::vector<double>> norm_log(norm.size(), inner);

	for (int i = 0; i != norm.size(); ++i) {
		for (int j = 0; j != norm[0].size(); ++j) {
			norm_log[i][j] = log(norm[i][j]);
		}
	}

	std::vector<double> slr_inner(2, 0.0);
	std::vector<std::vector<double>> slr(norm.size() - 1, slr_inner);
	for (int i = 0; i != norm.size() - 1; ++i) {
		slr[i] = regression::slr(norm_log[0], norm_log[i + 1]);
	}

	// Convergence rates.
	std::vector<double> 
		result{ slr[0][0], slr[1][0], slr[2][0], slr[3][0], slr[4][0] };

	if (show_output) {

		std::cout << std::scientific << std::setprecision(5);

		std::cout
			<< "SLR max-norm: " << std::endl
			<< std::setw(14) << slr[0][0]
			<< std::setw(14) << slr[0][1] << std::endl;

		std::cout <<
			"SLR l1 vector norm: " << std::endl
			<< std::setw(14) << slr[1][0]
			<< std::setw(14) << slr[1][1] << std::endl;

		std::cout <<
			"SLR l2 vector norm: " << std::endl
			<< std::setw(14) << slr[2][0]
			<< std::setw(14) << slr[2][1] << std::endl;

		std::cout <<
			"SLR L1 function norm: " << std::endl
			<< std::setw(14) << slr[3][0]
			<< std::setw(14) << slr[3][1] << std::endl;

		std::cout <<
			"SLR L2 function norm: " << std::endl
			<< std::setw(14) << slr[4][0]
			<< std::setw(14) << slr[4][1] << std::endl << std::endl;

		std::ofstream output_file(file_name);
		output_file << std::scientific << std::setprecision(12);
		for (int i = 0; i != norm[0].size(); ++i) {

			output_file
				<< std::setw(22) << norm[0][i] << ","
				<< std::setw(22) << norm_log[0][i] << ","
				<< std::setw(22) << norm[1][i] << ","
				<< std::setw(22) << norm_log[1][i] << ","
				<< std::setw(22) << norm[2][i] << ","
				<< std::setw(22) << norm_log[2][i] << ","
				<< std::setw(22) << norm[3][i] << ","
				<< std::setw(22) << norm_log[3][i] << ","
				<< std::setw(22) << norm[4][i] << ","
				<< std::setw(22) << norm_log[4][i] << ","
				<< std::setw(22) << norm[0][i] << ","
				<< std::setw(22) << norm_log[5][i] << std::endl;

		}
		output_file.close();

	}

	return result;

}


TEST(TriDiagonalSolver, HeatEquation1D) {

	// ####################################
	// Convergence rate in space dimension.
	// ####################################

	// Crank-Nicolson.
	{
		// Initial time grid.
		std::vector<double> time_grid = grid::exponential(0.0, 0.03, 101);

		// Initial spatial grid.
		std::vector<double> inner_spatial_grid = grid::hyperbolic(0.0, 1.0, 31);
		std::vector<std::vector<double>> spatial_grid(1, inner_spatial_grid);

		// Order of solution.
		const std::vector<int> inner_order(1, 1);
		const std::vector<std::vector<int>> order(1, inner_order);

		// Prefactors.
		const std::vector<double> inner_prefactor(1, 1.0);
		const std::vector<std::vector<double>> prefactor(1, inner_prefactor);

		// Diffusivity.
		const double diffusivity = 1.0;

		std::function<std::vector<double>
			(const double, const std::vector<std::vector<double>>&)>
			solution_generator = heat_eq::solution_func(
				order,
				prefactor,
				diffusivity
			);

		std::vector<std::vector<double>>
			norm = convergence::theta_1d<TriDiagonal>(
				time_grid,
				spatial_grid,
				grid::hyperbolic,
				d2dx2::nonuniform::c2b0,
				solution_generator,
				"space",
				20,
				5);

		std::vector<double> result = linear_regression(norm, true);

		// Maximum norm.
		EXPECT_NEAR(result[0], 2.0, 0.010);

		// L1 function norm.
		EXPECT_NEAR(result[3], 2.0, 0.013);

	}

	// Crank-Nicolson.
	{
		// Initial time grid.
		std::vector<double> time_grid = grid::uniform(0.0, 0.03, 201);

		// Initial spatial grid.
		std::vector<double> inner_spatial_grid = grid::hyperbolic(0.0, 1.0, 51);
		std::vector<std::vector<double>> spatial_grid(1, inner_spatial_grid);

		// Order of solution.
		const std::vector<int> inner_order(1, 3);
		const std::vector<std::vector<int>> order(1, inner_order);

		// Prefactors.
		const std::vector<double> inner_prefactor(1, 1.0);
		const std::vector<std::vector<double>> prefactor(1, inner_prefactor);

		// Diffusivity.
		const double diffusivity = 1.0;

		std::function<std::vector<double>
			(const double, const std::vector<std::vector<double>>&)>
			solution_generator = heat_eq::solution_func(
				order,
				prefactor,
				diffusivity
			);

		std::vector<std::vector<double>>
			norm = convergence::theta_1d<TriDiagonal>(
				time_grid,
				spatial_grid,
				grid::hyperbolic,
				d2dx2::nonuniform::c2b0,
				solution_generator,
				"space",
				20,
				5);

		std::vector<double> result = linear_regression(norm, true);

		// Maximum norm.
		EXPECT_NEAR(result[0], 2.0, 0.009);

		// L1 function norm.
		EXPECT_NEAR(result[3], 2.0, 0.005);

	}

	// Fully implicit.
	{
		// Initial time grid.
		std::vector<double> time_grid = grid::uniform(0.0, 0.03, 10001);

		// Initial spatial grid.
		std::vector<double> inner_spatial_grid = grid::hyperbolic(0.0, 1.0, 31);
		std::vector<std::vector<double>> spatial_grid(1, inner_spatial_grid);

		// Order of solution.
		const std::vector<int> inner_order(1, 1);
		const std::vector<std::vector<int>> order(1, inner_order);

		// Prefactors.
		const std::vector<double> inner_prefactor(1, 1.0);
		const std::vector<std::vector<double>> prefactor(1, inner_prefactor);

		// Diffusivity.
		const double diffusivity = 1.0;

		std::function<std::vector<double>
			(const double, const std::vector<std::vector<double>>&)>
			solution_generator = heat_eq::solution_func(
				order,
				prefactor,
				diffusivity
			);

		std::vector<std::vector<double>>
			norm = convergence::theta_1d<TriDiagonal>(
				time_grid,
				spatial_grid,
				grid::hyperbolic,
				d2dx2::nonuniform::c2b0,
				solution_generator,
				"space",
				20,
				5,
				1.0);

		std::vector<double> result = linear_regression(norm, true);

		// Maximum norm.
		EXPECT_NEAR(result[0], 2.0, 0.026);

		// L1 function norm.
		EXPECT_NEAR(result[3], 2.0, 0.036);

	}

	// ###################################
	// Convergence rate in time dimension.
	// ###################################

	// Crank-Nicolson.
	{
		// Initial time grid.
		std::vector<double> time_grid = grid::exponential(0.0, 0.03, 11);

		// Initial spatial grid.
		std::vector<double> inner_spatial_grid = grid::uniform(0.0, 1.0, 5001);
		std::vector<std::vector<double>> spatial_grid(1, inner_spatial_grid);

		// Order of solution.
		const std::vector<int> inner_order(1, 1);
		const std::vector<std::vector<int>> order(1, inner_order);

		// Prefactors.
		const std::vector<double> inner_prefactor(1, 1.0);
		const std::vector<std::vector<double>> prefactor(1, inner_prefactor);

		// Diffusivity.
		const double diffusivity = 1.0;

		std::function<std::vector<double>
			(const double, const std::vector<std::vector<double>>&)>
			solution_generator = heat_eq::solution_func(
				order,
				prefactor,
				diffusivity
			);

		std::vector<std::vector<double>>
			norm = convergence::theta_1d<TriDiagonal>(
				time_grid,
				spatial_grid,
				grid::exponential,
				d2dx2::uniform::c2b0,
				solution_generator,
				"time",
				20,
				5);

		std::vector<double> result = linear_regression(norm, true);

		// Maximum norm.
		EXPECT_NEAR(result[0], 2.0, 0.008);

		// L1 function norm.
		EXPECT_NEAR(result[3], 2.0, 0.008);

	}

	// Crank-Nicolson.
	{
		// Initial time grid.
		std::vector<double> time_grid = grid::uniform(0.0, 0.03, 11);

		// Initial spatial grid.
		std::vector<double> inner_spatial_grid = grid::uniform(0.0, 1.0, 5001);
		std::vector<std::vector<double>> spatial_grid(1, inner_spatial_grid);

		// Order of solution.
		const std::vector<int> inner_order(1, 3);
		const std::vector<std::vector<int>> order(1, inner_order);

		// Prefactors.
		const std::vector<double> inner_prefactor(1, 1.0);
		const std::vector<std::vector<double>> prefactor(1, inner_prefactor);

		// Diffusivity.
		const double diffusivity = 1.0;

		std::function<std::vector<double>
			(const double, const std::vector<std::vector<double>>&)>
			solution_generator = heat_eq::solution_func(
				order,
				prefactor,
				diffusivity
			);

		std::vector<std::vector<double>>
			norm = convergence::theta_1d<TriDiagonal>(
				time_grid,
				spatial_grid,
				grid::uniform,
				d2dx2::uniform::c2b0,
				solution_generator,
				"time",
				20,
				5);

		std::vector<double> result = linear_regression(norm, true);

		// Maximum norm.
		EXPECT_NEAR(result[0], 2.0, 0.004);

		// L1 function norm.
		EXPECT_NEAR(result[3], 2.0, 0.004);

	}

	// Fully implicit.
	{
		// Initial time grid.
		std::vector<double> time_grid = grid::uniform(0.0, 0.03, 11);

		// Initial spatial grid.
		std::vector<double> inner_spatial_grid = grid::uniform(0.0, 1.0, 5001);
		std::vector<std::vector<double>> spatial_grid(1, inner_spatial_grid);

		// Order of solution.
		const std::vector<int> inner_order(1, 1);
		const std::vector<std::vector<int>> order(1, inner_order);

		// Prefactors.
		const std::vector<double> inner_prefactor(1, 1.0);
		const std::vector<std::vector<double>> prefactor(1, inner_prefactor);

		// Diffusivity.
		const double diffusivity = 1.0;

		std::function<std::vector<double>
			(const double, const std::vector<std::vector<double>>&)>
			solution_generator = heat_eq::solution_func(
				order,
				prefactor,
				diffusivity
			);

		std::vector<std::vector<double>>
			norm = convergence::theta_1d<TriDiagonal>(
				time_grid,
				spatial_grid,
				grid::uniform,
				d2dx2::uniform::c2b0,
				solution_generator,
				"time",
				20,
				5,
				1.0);

		std::vector<double> result = linear_regression(norm, true);

		// Maximum norm.
		EXPECT_NEAR(result[0], 1.0, 0.057);

		// L1 function norm.
		EXPECT_NEAR(result[3], 1.0, 0.057);

	}

}


TEST(PentaDiagonalSolver, HeatEquation1D) {

	// Crank-Nicolson.
	{
		// Initial time grid.
		std::vector<double> time_grid = grid::uniform(0.0, 0.03, 501);

		// Initial spatial grid.
		std::vector<double> inner_spatial_grid = grid::uniform(0.0, 1.0, 11);
		std::vector<std::vector<double>> spatial_grid(1, inner_spatial_grid);

		// Order of solution.
		const std::vector<int> inner_order(1, 1);
		const std::vector<std::vector<int>> order(1, inner_order);

		// Prefactors.
		const std::vector<double> inner_prefactor(1, 1.0);
		const std::vector<std::vector<double>> prefactor(1, inner_prefactor);

		// Diffusivity.
		const double diffusivity = 1.0;

		std::function<std::vector<double>
			(const double, const std::vector<std::vector<double>>&)>
			solution_generator = heat_eq::solution_func(
				order,
				prefactor,
				diffusivity
			);

		std::vector<std::vector<double>>
			norm = convergence::theta_1d<PentaDiagonal>(
				time_grid,
				spatial_grid,
				grid::uniform,
				d2dx2::uniform::c4b0,
				solution_generator,
				"space",
				11,
				5);

		std::vector<double> result = linear_regression(norm, true);

		// Maximum norm.
		EXPECT_NEAR(result[0], 4.8, 0.017);

		// L1 function norm.
		EXPECT_NEAR(result[3], 4.6, 0.059);

	}

	// Crank-Nicolson.
	{
		// Initial time grid.
		std::vector<double> time_grid = grid::uniform(0.0, 0.03, 1001);

		// Initial spatial grid.
		std::vector<double> inner_spatial_grid = grid::hyperbolic(0.0, 1.0, 21);
		std::vector<std::vector<double>> spatial_grid(1, inner_spatial_grid);

		// Order of solution.
		const std::vector<int> inner_order(1, 1);
		const std::vector<std::vector<int>> order(1, inner_order);

		// Prefactors.
		const std::vector<double> inner_prefactor(1, 1.0);
		const std::vector<std::vector<double>> prefactor(1, inner_prefactor);

		// Diffusivity.
		const double diffusivity = 1.0;

		std::function<std::vector<double>
			(const double, const std::vector<std::vector<double>>&)>
			solution_generator = heat_eq::solution_func(
				order,
				prefactor,
				diffusivity
			);

		std::vector<std::vector<double>>
			norm = convergence::theta_1d<PentaDiagonal>(
				time_grid,
				spatial_grid,
				grid::hyperbolic,
				d2dx2::nonuniform::c4b0,
				solution_generator,
				"space",
				11,
				5);

		std::vector<double> result = linear_regression(norm, true);

		// Maximum norm.
		EXPECT_NEAR(result[0], 3.5, 0.065);

		// L1 function norm.
		EXPECT_NEAR(result[3], 3.6, 0.009);

	}

}


TEST(TriDiagonalSolver, HeatEquation2D) {

	// Crank-Nicolson.
	{
		// Initial time grid.
		std::vector<double> time_grid = grid::uniform(0.0, 0.03, 201);

		// Initial spatial grid.
		std::vector<double> spatial_grid_x = grid::uniform(0.0, 1.0, 11);
		std::vector<double> spatial_grid_y = grid::uniform(0.0, 1.0, 201);
		std::vector<std::vector<double>> spatial_grid {spatial_grid_x, spatial_grid_y};

		// Order of solution.
		const std::vector<int> inner_order(2, 1);
		const std::vector<std::vector<int>> order(1, inner_order);

		// Prefactors.
		const std::vector<double> inner_prefactor(2, 1.0);
		const std::vector<std::vector<double>> prefactor(1, inner_prefactor);

		// Diffusivity.
		const double diffusivity = 1.0;

		std::function<std::vector<double>
			(const double, const std::vector<std::vector<double>>&)>
			solution_generator = heat_eq::solution_func(
				order,
				prefactor,
				diffusivity
			);

		std::vector<std::vector<double>>
			norm = convergence::adi::dr_2d<TriDiagonal, TriDiagonal>(
				time_grid,
				spatial_grid,
				grid::uniform,
				d2dx2::uniform::c2b0,
				d2dx2::uniform::c2b0,
				solution_generator,
				"space_1",
				11,
				5);

		std::vector<double> result = linear_regression(norm, true);

		// Maximum norm.
		EXPECT_NEAR(result[0], 2.0, 0.047);

		// L1 function norm.
		EXPECT_NEAR(result[3], 2.0, 0.050);

	}

	// Crank-Nicolson.
	{
		// Initial time grid.
		std::vector<double> time_grid = grid::uniform(0.0, 0.03, 201);

		// Initial spatial grid.
		std::vector<double> spatial_grid_x = grid::uniform(0.0, 1.0, 201);
		std::vector<double> spatial_grid_y = grid::uniform(0.0, 1.0, 11);
		std::vector<std::vector<double>> spatial_grid{ spatial_grid_x, spatial_grid_y };

		// Order of solution.
		const std::vector<int> inner_order(2, 1);
		const std::vector<std::vector<int>> order(1, inner_order);

		// Prefactors.
		const std::vector<double> inner_prefactor(2, 1.0);
		const std::vector<std::vector<double>> prefactor(1, inner_prefactor);

		// Diffusivity.
		const double diffusivity = 1.0;

		std::function<std::vector<double>
			(const double, const std::vector<std::vector<double>>&)>
			solution_generator = heat_eq::solution_func(
				order,
				prefactor,
				diffusivity
			);

		std::vector<std::vector<double>>
			norm = convergence::adi::dr_2d<TriDiagonal, TriDiagonal>(
				time_grid,
				spatial_grid,
				grid::uniform,
				d2dx2::uniform::c2b0,
				d2dx2::uniform::c2b0,
				solution_generator,
				"space_2",
				11,
				5);

		std::vector<double> result = linear_regression(norm, true);

		// Maximum norm.
		EXPECT_NEAR(result[0], 2.0, 0.047);

		// L1 function norm.
		EXPECT_NEAR(result[3], 2.0, 0.060);

	}

	// Crank-Nicolson.
	{
		// Initial time grid.
		std::vector<double> time_grid = grid::uniform(0.0, 0.03, 21);

		// Initial spatial grid.
		std::vector<double> spatial_grid_x = grid::uniform(0.0, 1.0, 1501);
		std::vector<double> spatial_grid_y = grid::uniform(0.0, 1.0, 1501);
		std::vector<std::vector<double>> spatial_grid{ spatial_grid_x, spatial_grid_y };

		// Order of solution.
		const std::vector<int> inner_order(2, 1);
		const std::vector<std::vector<int>> order(1, inner_order);

		// Prefactors.
		const std::vector<double> inner_prefactor(2, 1.0);
		const std::vector<std::vector<double>> prefactor(1, inner_prefactor);

		// Diffusivity.
		const double diffusivity = 1.0;

		std::function<std::vector<double>
			(const double, const std::vector<std::vector<double>>&)>
			solution_generator = heat_eq::solution_func(
				order,
				prefactor,
				diffusivity
			);

		std::vector<std::vector<double>>
			norm = convergence::adi::dr_2d<TriDiagonal, TriDiagonal>(
				time_grid,
				spatial_grid,
				grid::uniform,
				d2dx2::uniform::c2b0,
				d2dx2::uniform::c2b0,
				solution_generator,
				"time",
				3,
				5);

		std::vector<double> result = linear_regression(norm, true);

		// Maximum norm.
		EXPECT_NEAR(result[0], 2.0, 0.064);

		// L1 function norm.
		EXPECT_NEAR(result[3], 2.0, 0.064);

	}

}


TEST(TriDiagonalSolver, HeatEquation2D_new) {

	// Crank-Nicolson.
	{
		// Initial time grid.
		std::vector<double> time_grid = grid::uniform(0.0, 0.03, 201);

		// Initial spatial grid.
		std::vector<double> spatial_grid_x = grid::uniform(0.0, 1.0, 11);
		std::vector<double> spatial_grid_y = grid::uniform(0.0, 1.0, 201);
		std::vector<std::vector<double>> spatial_grid{ spatial_grid_x, spatial_grid_y };

		// Order of solution.
		const std::vector<int> inner_order(2, 1);
		const std::vector<std::vector<int>> order(1, inner_order);

		// Prefactors.
		const std::vector<double> inner_prefactor(2, 1.0);
		const std::vector<std::vector<double>> prefactor(1, inner_prefactor);

		// Diffusivity.
		const double diffusivity = 1.0;

		std::function<std::vector<double>
			(const double, const std::vector<std::vector<double>>&)>
			solution_generator = heat_eq::solution_func(
				order,
				prefactor,
				diffusivity
			);

		// SHOULD CALL OVERLOADED dr_2d with prefactors and vectors of derivative_generators!!!

		std::vector<std::vector<double>>
			norm = convergence::adi::dr_2d<TriDiagonal, TriDiagonal>(
				time_grid,
				spatial_grid,
				grid::uniform,
				d2dx2::uniform::c2b0,
				d2dx2::uniform::c2b0,
				solution_generator,
				"space_1",
				11,
				5);

		std::vector<double> result = linear_regression(norm, true);

		// Maximum norm.
		EXPECT_NEAR(result[0], 2.0, 0.047);

		// L1 function norm.
		EXPECT_NEAR(result[3], 2.0, 0.050);

	}

}
