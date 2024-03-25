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


std::vector<std::vector<double>> prefactor_generator(
	const std::vector<std::vector<double>>& grid);


// Prefactors C * f(grid_1) * g(grid_2) * ... for 1st and 2nd order derivatives.
std::vector<std::vector<double>> prefactor_generator(
	const std::vector<std::vector<double>>& grid) {

	std::vector<std::vector<double>> result(3, { 1.0 });

	// #######################
	// First order derivative.
	// #######################
	
	result[0][0] = 0.0;

	for (int i = 0; i != grid.size(); ++i) {
		for (int j = 0; j != grid[i].size(); ++j) {

			result[0].push_back(0.0);
//			result[0].push_back(grid[i][j]);

		}
	}
			
	// ########################
	// Second order derivative.
	// ########################

	result[1][0] = 1.0;

	for (int i = 0; i != grid.size(); ++i) {
		for (int j = 0; j != grid[i].size(); ++j) {

			result[1].push_back(1.0);
//			result[1].push_back(grid[i][j]);

		}
	}

	// ###################
	// Inhomogenious term.
	// ###################

	result[2][0] = 0.0;

	for (int i = 0; i != grid.size(); ++i) {
		for (int j = 0; j != grid[i].size(); ++j) {
			result[2].push_back(0.0);
		}
	}

	return result;

}


std::vector<double> mixed_prefactor_generator(
	const std::vector<std::vector<double>>& grid);


// Prefactors C * f(grid_1) * g(grid_2) for mixed derivative.
std::vector<double> mixed_prefactor_generator(
	const std::vector<std::vector<double>>& grid) {

	int n_points = grid[0].size() * grid[1].size();

	std::vector<double> result(n_points, 0.0);

	int index = 0;
	for (int i = 0; i != grid[0].size(); ++i) {
		for (int j = 0; j != grid[1].size(); ++j) {
			result[index] = 0.0;
			++index;
		}
	}

	return result;

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


		std::vector<std::function<TriDiagonal(std::vector<double>)>> 
			deriv_1{ d1dx1::uniform::c2b1, d2dx2::uniform::c2b0 };
		std::vector<std::function<TriDiagonal(std::vector<double>)>>
			deriv_2{ d1dx1::uniform::c2b1, d2dx2::uniform::c2b0 };

#if false
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
#endif
#if false
		std::vector<std::vector<double>>
			norm = convergence::adi::dr_2d<TriDiagonal, TriDiagonal>(
				time_grid,
				spatial_grid,
				grid::uniform,
				{ 0.0, diffusivity },
				{ 0.0, diffusivity },
				deriv_1,
				deriv_2,
				solution_generator,
				"space_1",
				11,
				5);
#endif
#if false
		std::vector<std::vector<double>>
			norm = convergence::adi::dr_2d<TriDiagonal, TriDiagonal>(
				time_grid,
				spatial_grid,
				grid::uniform,
				prefactor_generator,
				prefactor_generator,
				deriv_1,
				deriv_2,
				solution_generator,
				"space_1",
				11,
				5);
#endif
#if false
		std::vector<std::vector<double>>
			norm = convergence::adi::cs_2d<TriDiagonal, TriDiagonal>(
				time_grid,
				spatial_grid,
				grid::uniform,
				d2dx2::uniform::c2b0,
				d2dx2::uniform::c2b0,
				solution_generator,
				"space_1",
				11,
				5);
#endif
#if false
		std::vector<std::vector<double>>
			norm = convergence::adi::cs_2d<TriDiagonal, TriDiagonal>(
				time_grid,
				spatial_grid,
				grid::uniform,
				{ 0.0, diffusivity },
				{ 0.0, diffusivity },
				deriv_1,
				deriv_2,
				solution_generator,
				"space_1",
				11,
				5);
#endif
		
		std::vector<std::vector<double>>
			norm = convergence::adi::cs_2d<TriDiagonal, TriDiagonal>(
				time_grid,
				spatial_grid,
				grid::uniform,
				prefactor_generator,
				prefactor_generator,
				mixed_prefactor_generator,
				deriv_1,
				deriv_2,
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



// Prefactors C * f(grid_1) * g(grid_2) * ... for 1st and 2nd order derivatives.
std::vector<std::vector<double>> prefactor_generator_heston_s(
	const std::vector<std::vector<double>>& grid,
	const double rate) {

	std::vector<std::vector<double>> result(3, { 1.0 });

	// #######################
	// First order derivative.
	// #######################

	result[0][0] = rate;

	for (int i = 0; i != grid.size(); ++i) {
		for (int j = 0; j != grid[i].size(); ++j) {
			if (i == 0) {
				result[0].push_back(grid[i][j]);
			}
			else {
				result[0].push_back(1.0);
			}
		}
	}

	// ########################
	// Second order derivative.
	// ########################

	result[1][0] = 0.5;

	for (int i = 0; i != grid.size(); ++i) {
		for (int j = 0; j != grid[i].size(); ++j) {
			if (i == 0) {
				result[1].push_back(grid[i][j] * grid[i][j]);
			}
			else {
				result[1].push_back(grid[i][j]);
			}
		}
	}

	// ###################
	// Inhomogenious term.
	// ###################

	result[2][0] = -0.5 * rate;

	for (int i = 0; i != grid.size(); ++i) {
		for (int j = 0; j != grid[i].size(); ++j) {
			result[2].push_back(1.0);
		}
	}

	return result;

}

// Prefactors C * f(grid_1) * g(grid_2) * ... for 1st and 2nd order derivatives.
std::vector<std::vector<double>> prefactor_generator_heston_v(
	const std::vector<std::vector<double>>& grid,
	const double rate,
	const double lambda,
	const double theta,
	const double eta) {

	std::vector<std::vector<double>> result(3, { 1.0 });

	// #######################
	// First order derivative.
	// #######################

	result[0][0] = lambda;

	for (int i = 0; i != grid.size(); ++i) {
		for (int j = 0; j != grid[i].size(); ++j) {
			if (i == 1) {
				result[0].push_back(theta - grid[i][j]);
			}
			else {
				result[0].push_back(1.0);
			}
		}
	}

	// ########################
	// Second order derivative.
	// ########################

	result[1][0] = 0.5 * eta * eta;

	for (int i = 0; i != grid.size(); ++i) {
		for (int j = 0; j != grid[i].size(); ++j) {
			if (i == 1) {
				result[1].push_back(grid[i][j]);
			}
			else {
				result[1].push_back(1.0);
			}
		}
	}

	// ###################
	// Inhomogenious term.
	// ###################

	result[2][0] = -0.5 * rate;

	for (int i = 0; i != grid.size(); ++i) {
		for (int j = 0; j != grid[i].size(); ++j) {
			result[2].push_back(1.0);
		}
	}

	return result;

}

// Prefactors C * f(grid_1) * g(grid_2) for mixed derivative.
std::vector<double> mixed_prefactor_generator_heston(
	const std::vector<std::vector<double>>& grid,
	const double eta,
	const double rho) {

	int n_points = 1;
	for (int i = 0; i != grid.size(); ++i) {
		n_points *= grid[i].size();
	}

	std::vector<double> result(n_points, 0.0);

	int index = 0;
	for (int i = 0; i != grid[0].size(); ++i) {
		for (int j = 0; j != grid[1].size(); ++j) {
			result[index] = eta * rho * grid[0][i] * grid[1][j];
			++index;
		}
	}

	return result;

}


TEST(TriDiagonalSolver, HestonCall) {

	const double rate = 0.03;

	const double lambda = 3.0;
	const double theta = 0.12;
	const double eta = 0.041;

	const double rho = 0.6;

	const double strike = 100.0;

	// Crank-Nicolson.
	{

		// Initial time grid.
		std::vector<double> time_grid = grid::uniform(0.0, 1.0, 31);

		// Initial spatial grid.
		std::vector<double> spatial_grid_x = grid::uniform(2.0, 200.0, 51);
		std::vector<double> spatial_grid_y = grid::uniform(0.0, 1.0, 21);
		std::vector<std::vector<double>> spatial_grid{ spatial_grid_x, spatial_grid_y };


		std::vector<std::vector<double>> prefactors_1 = prefactor_generator_heston_s(spatial_grid, rate);
		std::vector<std::vector<double>> prefactors_2 = prefactor_generator_heston_v(spatial_grid, rate, lambda, theta, eta);


		std::vector<std::function<TriDiagonal(std::vector<double>)>>
			deriv_1{ d1dx1::uniform::c2b1, d2dx2::uniform::c2b0 };
		std::vector<std::function<TriDiagonal(std::vector<double>)>>
			deriv_2{ d1dx1::uniform::c2b1, d2dx2::uniform::c2b0 };

		std::vector<TriDiagonal> derivatives_1;
		derivatives_1.push_back(deriv_1[0](spatial_grid[0]).identity());
		derivatives_1.push_back(deriv_1[0](spatial_grid[0]));
		derivatives_1.push_back(deriv_1[1](spatial_grid[0]));

		std::vector<TriDiagonal> derivatives_2;
		derivatives_2.push_back(deriv_2[0](spatial_grid[1]).identity());
		derivatives_2.push_back(deriv_2[0](spatial_grid[1]));
		derivatives_2.push_back(deriv_2[1](spatial_grid[1]));


		std::vector<double> prefactors_12 = mixed_prefactor_generator_heston(spatial_grid, eta, rho);


		MixedDerivative<TriDiagonal, TriDiagonal> mixed(
			deriv_1[0](spatial_grid[0]),
			deriv_2[0](spatial_grid[1]));


		mixed.set_prefactors(prefactors_12);


		std::vector<double> func(spatial_grid_x.size() * spatial_grid_y.size(), 0.0);
		int index = 0;
		for (int i = 0; i != spatial_grid_x.size(); ++i) {
			for (int j = 0; j != spatial_grid_y.size(); ++j) {
				func[index] = std::max(spatial_grid_x[i] - strike, 0.0);
				++index;
			}
		}


		propagation::adi::cs_2d(
			time_grid,
			prefactors_1, 
			prefactors_2,
			derivatives_1, 
			derivatives_2,
			mixed,
			func,
			0.5);


		std::ofstream file("heston_surface.csv");
		file << std::scientific << std::setprecision(12);
		index = 0;
		for (int i = 0; i != spatial_grid_x.size(); ++i) {
			for (int j = 0; j != spatial_grid_y.size(); ++j) {
				file << std::setw(22) << func[index];
				if (j < spatial_grid_y.size() - 1) {
					file << ", ";
				}
				++index;
			}
			file << std::endl;
		}
		file << std::endl;
		file.close();

		std::ofstream file_s("heston_surface_s.csv");
		file_s << std::scientific << std::setprecision(12);
		index = 0;
		for (int i = 0; i != spatial_grid_x.size(); ++i) {
			file_s << std::setw(22) << spatial_grid_x[i];
			if (i < spatial_grid_x.size() - 1) {
				file_s << ", ";
			}
		}
		file_s << std::endl;
		file_s.close();

		std::ofstream file_v("heston_surface_v.csv");
		file_v << std::scientific << std::setprecision(12);
		index = 0;
		for (int i = 0; i != spatial_grid_y.size(); ++i) {
			file_v << std::setw(22) << spatial_grid_y[i];
			if (i < spatial_grid_y.size() - 1) {
				file_v << ", ";
			}
		}
		file_v << std::endl;
		file_v.close();

	}

}


TEST(TriDiagonalSolver, BlackScholesCall) {

	const double rate = 0.03;

	const double sigma = 0.2;

	const double tau = 1.0;

	const double strike = 100.0;

	// Crank-Nicolson.
	{

		// Initial time grid.
		std::vector<double> time_grid = grid::uniform(0.0, tau, 31);

		// Initial spatial grid.
		std::vector<double> spatial_grid = grid::uniform(0.0, 200.0, 51);

		std::vector<double> prefactor_1(spatial_grid.size(), 0.0);
		std::vector<double> prefactor_2(spatial_grid.size(), 0.0);

		for (int i = 0; i != spatial_grid.size(); ++i) {
			prefactor_1[i] = rate * spatial_grid[i];
			prefactor_2[i] = 0.5 * sigma * sigma * spatial_grid[i] * spatial_grid[i];
		}

		std::vector<std::function<TriDiagonal(std::vector<double>)>>
			deriv{ d1dx1::uniform::c2b1, d2dx2::uniform::c2b0 };

		TriDiagonal test_1 = deriv[0](spatial_grid);
		TriDiagonal test_2 = deriv[0](spatial_grid).pre_vector(prefactor_1);
		TriDiagonal test_3 = deriv[1](spatial_grid);
		TriDiagonal test_4 = deriv[1](spatial_grid).pre_vector(prefactor_2);

		TriDiagonal derivative = deriv[0](spatial_grid).pre_vector(prefactor_1);
		derivative += deriv[1](spatial_grid).pre_vector(prefactor_2);
		derivative -= rate * deriv[0](spatial_grid).identity();

		std::vector<double> func(spatial_grid.size(), 0.0);
		for (int i = 0; i != spatial_grid.size(); ++i) {
			func[i] = std::max(spatial_grid[i] - strike, 0.0);
		}
		std::vector<double> ic_func = func;

		propagation::theta_1d::full(
			time_grid,
			derivative,
			func,
			0.5);

		std::ofstream file("black_scholes_call.csv");
		file << std::scientific << std::setprecision(12);
		for (int i = 0; i != spatial_grid.size(); ++i) {
			file 
				<< std::setw(22) << spatial_grid[i] << ","
				<< std::setw(22) << func[i] << ","
				<< std::setw(22) << bs::call::price(spatial_grid[i], rate, sigma, strike, tau) << ","
				<< std::setw(22) << ic_func[i]
				<< std::endl;
		}
		file.close();

	}

}
