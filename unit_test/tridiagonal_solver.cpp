#include "pch.h"


std::vector<double> linear_regression(
	std::vector<double>& step_size_vec,
	std::vector<double>& max_norm,
	std::vector<double>& l1_vec_norm,
	std::vector<double>& l2_vec_norm,
	std::vector<double>& l1_func_norm,
	std::vector<double>& l2_func_norm,
	const bool show_output = true,
	const std::string file_name = "slr.csv") {

	std::vector<double> step_size_vec_log(step_size_vec.size(), 0.0);
	std::vector<double> max_norm_log(max_norm.size(), 0.0);
	std::vector<double> l1_vec_norm_log(l1_vec_norm.size(), 0.0);
	std::vector<double> l2_vec_norm_log(l2_vec_norm.size(), 0.0);
	std::vector<double> l1_func_norm_log(l1_func_norm.size(), 0.0);
	std::vector<double> l2_func_norm_log(l2_func_norm.size(), 0.0);

	for (int i = 0; i != step_size_vec.size(); ++i) {

		step_size_vec_log[i] = log(step_size_vec[i]);
		max_norm_log[i] = log(max_norm[i]);
		l1_vec_norm_log[i] = log(l1_vec_norm[i]);
		l2_vec_norm_log[i] = log(l2_vec_norm[i]);
		l1_func_norm_log[i] = log(l1_func_norm[i]);
		l2_func_norm_log[i] = log(l2_func_norm[i]);

	}

	std::vector<double> slr_max = regression::slr(step_size_vec_log, max_norm_log);
	std::vector<double> slr_l1_vec = regression::slr(step_size_vec_log, l1_vec_norm_log);
	std::vector<double> slr_l2_vec = regression::slr(step_size_vec_log, l2_vec_norm_log);
	std::vector<double> slr_l1_func = regression::slr(step_size_vec_log, l1_func_norm_log);
	std::vector<double> slr_l2_func = regression::slr(step_size_vec_log, l2_func_norm_log);

	// Convergence rates.
	std::vector<double> result{ slr_max[0], slr_l1_vec[0], slr_l2_vec[0], slr_l1_func[0], slr_l2_func[0] };

	if (show_output) {

		std::cout << std::scientific << std::setprecision(5);

		std::cout
			<< "SLR max-norm: " << std::endl
			<< std::setw(14) << slr_max[0]
			<< std::setw(14) << slr_max[1] << std::endl;

		std::cout <<
			"SLR l1 vector norm: " << std::endl
			<< std::setw(14) << slr_l1_vec[0]
			<< std::setw(14) << slr_l1_vec[1] << std::endl;

		std::cout <<
			"SLR l2 vector norm: " << std::endl
			<< std::setw(14) << slr_l2_vec[0]
			<< std::setw(14) << slr_l2_vec[1] << std::endl;

		std::cout <<
			"SLR L1 function norm: " << std::endl
			<< std::setw(14) << slr_l1_func[0]
			<< std::setw(14) << slr_l1_func[1] << std::endl;

		std::cout <<
			"SLR L2 function norm: " << std::endl
			<< std::setw(14) << slr_l2_func[0]
			<< std::setw(14) << slr_l2_func[1] << std::endl << std::endl;

		std::ofstream output_file(file_name);
		output_file << std::scientific << std::setprecision(12);
		for (int i = 0; i != step_size_vec.size(); ++i) {

			output_file
				<< std::setw(22) << step_size_vec[i] << ","
				<< std::setw(22) << step_size_vec_log[i] << ","
				<< std::setw(22) << max_norm[i] << ","
				<< std::setw(22) << max_norm_log[i] << ","
				<< std::setw(22) << l1_vec_norm[i] << ","
				<< std::setw(22) << l1_vec_norm_log[i] << ","
				<< std::setw(22) << l2_vec_norm[i] << ","
				<< std::setw(22) << l2_vec_norm_log[i] << ","
				<< std::setw(22) << l1_func_norm[i] << ","
				<< std::setw(22) << l1_func_norm_log[i] << ","
				<< std::setw(22) << l2_func_norm[i] << ","
				<< std::setw(22) << l2_func_norm_log[i] << std::endl;

		}
		output_file.close();

	}

	return result;

}


std::vector<double> heat_equation_1d(
	const std::string convergence_type,
	const int x_points_start,
	const std::string x_grid_type,
	const int t_points_start,
	const double t_interval,
	const std::string t_grid_type,
	const std::string d2dx2_type,
	const std::string solution_type,
	const int n_iterations = 50,
	const double theta = 0.5) {

	bool show_output_basic = true;
	bool show_output_all = false;

	const int n_increments = 5;

	std::vector<double> step_size_vec;
	std::vector<double> max_norm;
	std::vector<double> l1_vec_norm;
	std::vector<double> l2_vec_norm;
	std::vector<double> l1_func_norm;
	std::vector<double> l2_func_norm;

	int t_points = t_points_start;
	int x_points = x_points_start;

	// Lambda parameter. See wiki page for heat equation.
	double lambda = 0.0;
	if (solution_type == "cos(pi*x)") {
		lambda = pow(1 * M_PI, 2);
	}
	else if (solution_type == "cos(3*pi*x)") {
		lambda = pow(3 * M_PI, 2);
	}
	else {
		throw std::invalid_argument("Unknown solution type.");
	}

	// Convergence rate wrt spatial dimension.
	for (int i = 0; i != n_iterations; ++i) {

		if (convergence_type == "time") {
			// Number of grid points in time dimension.
			t_points = t_points_start + n_increments * i;
		}
		else if (convergence_type == "space") {
			// Number of grid points in spatial dimension.
			x_points = x_points_start + n_increments * i;
		}
		else {
			throw std::invalid_argument("convergence_type unknown.");
		}

		// Time grid.
		std::vector<double> t_grid;
		if (t_grid_type == "uniform")
		{
			t_grid = grid::uniform(0.0, t_interval, t_points);
		}
		else if (t_grid_type == "exponential")
		{
			t_grid = grid::exponential(0.0, t_interval, t_points);
		}
		else if (t_grid_type == "hyperbolic")
		{
			t_grid = grid::hyperbolic(0.0, t_interval, t_points);
		}
		else {
			throw std::invalid_argument("t_grid_type unknown.");
		}

		// Average time grid spacing.
		const double dt = t_interval / (t_points - 1);

		// Spatial grid. TODO: Chosen such that d2dx2 = 0 at boundary for analytical solution!
		std::vector<double> x_grid;
		if (x_grid_type == "uniform")
		{
			x_grid = grid::uniform(-0.5, 0.5, x_points);
		}
		else if (x_grid_type == "exponential")
		{
			x_grid = grid::exponential(-0.5, 0.5, x_points);
		}
		else if (x_grid_type == "hyperbolic")
		{
			x_grid = grid::hyperbolic(-0.5, 0.5, x_points);
		}
		else {
			throw std::invalid_argument("x_grid_type unknown.");
		}

		// Average spatial grid spacing.
		const double dx = (x_grid.back() - x_grid.front()) / (x_points - 1);

		if (convergence_type == "time") {
			step_size_vec.push_back(dt);
		}
		else if (convergence_type == "space") {
			step_size_vec.push_back(dx);
		}
		else {
			throw std::invalid_argument("convergence_type unknown.");
		}

		// Initial condition.
		std::vector<double> func(x_points, 0.0);
		for (int i = 0; i != x_points; ++i) {
			func[i] = cos(sqrt(lambda) * x_grid[i]);
		}

		// ################
		// Setup operators.
		// ################

		// Second order differential operator.
		TriDiagonal deriv_operator(x_points);
		PentaDiagonal deriv_operator_p(x_points);

		if (d2dx2_type == "d2dx2::uniform::c2b0") {
			deriv_operator = d2dx2::uniform::c2b0(x_points, dx);
		}

		else if (d2dx2_type == "d2dx2::uniform::c4b0") {
			deriv_operator_p = d2dx2::uniform::c4b0(x_points, dx);
		}

		else if (d2dx2_type == "d2dx2::nonuniform::c2b0") {
			deriv_operator = d2dx2::nonuniform::c2b0(x_points, x_grid);
		}

		else if (d2dx2_type == "d2dx2::nonuniform::c4b0") {
			deriv_operator_p = d2dx2::nonuniform::c4b0(x_points, x_grid);
		}

		else {
			throw std::invalid_argument("d2dx2_type unknown.");
		}

		// Identity operator.
//		TriDiagonal iden = identity::tri(x_points, deriv_operator.n_boundary_elements());
		TriDiagonal iden = deriv_operator.identity();
//		PentaDiagonal iden_p = identity::penta(x_points, deriv_operator_p.n_boundary_elements() - 1);
		PentaDiagonal iden_p = deriv_operator_p.identity();

		// ############
		// FD solution.
		// ############

		std::vector<double> solution = func;

		if (d2dx2_type != "d2dx2::uniform::c4b0" && d2dx2_type != "d2dx2::nonuniform::c4b0") {

			TriDiagonal lhs(x_points);
			TriDiagonal rhs(x_points);

			for (int i = 0; i != t_points - 1; ++i) {

				double dt_tmp = t_grid[i + 1] - t_grid[i];

				// LHS operator.
				lhs = deriv_operator;
				lhs *= -theta * dt_tmp;
				lhs += iden;

				// RHS operator.
				rhs = deriv_operator;
				rhs *= (1.0 - theta) * dt_tmp;
				rhs += iden;

				// Evaluation RHS.
				solution = rhs * solution;

				// Solve matrix equation.
				tri_solver(lhs, solution);

			}
		}
		else {

			PentaDiagonal lhs(x_points, 2, deriv_operator_p.n_boundary_elements() - 1);
			PentaDiagonal rhs(x_points, 2, deriv_operator_p.n_boundary_elements() - 1);

			for (int i = 0; i != t_points - 1; ++i) {

				double dt_tmp = t_grid[i + 1] - t_grid[i];

				// LHS operator.
				lhs = deriv_operator_p;
				lhs *= -theta * dt_tmp;
				lhs += iden_p;

				// RHS operator.
				rhs = deriv_operator_p;
				rhs *= (1.0 - theta) * dt_tmp;
				rhs += iden_p;

				// Evaluation RHS.
				solution = rhs * solution;

				// Solve matrix equation.
				penta_solver(lhs, solution);

			}

		}

		// ####################
		// Analytical solution.
		// ####################

		std::vector<double> analytical(x_points, 0.0);
		for (int i = 0; i != x_points; ++i) {
			analytical[i] = exp(-lambda * t_interval) * func[i];
		}

		if (show_output_all) {
			std::cout << std::scientific << std::setprecision(5);
			std::cout << "Average time step: " << dt << "\t Time interval: " << (t_grid.back() - t_grid.front()) << std::endl;
			for (int i = 0; i != x_points; ++i) {
				std::cout
					<< std::setw(3) << i
					<< std::setw(14) << x_grid[i]
					<< std::setw(14) << func[i]
					<< std::setw(14) << analytical[i]
					<< std::setw(14) << solution[i]
					<< std::setw(14) << abs(analytical[i] - solution[i])
					<< std::endl;

			}
			std::cout << std::endl;
		}

		// Difference vector.
		std::vector<double> diff = norm::vector_diff(analytical, solution);

		max_norm.push_back(norm::vector::max(diff));

		l1_vec_norm.push_back(norm::vector::l1(diff));

		l2_vec_norm.push_back(norm::vector::l2(diff));

		l1_func_norm.push_back(norm::function::l1(x_grid, diff));

		l2_func_norm.push_back(norm::function::l2(x_grid, diff));

	}

	std::vector<double> result = linear_regression(step_size_vec, max_norm, l1_vec_norm, l2_vec_norm, l1_func_norm, l2_func_norm, true);

	return result;

}


TEST(TriDiagonalSolver, HeatEquation1D) {

	// ####################################
	// Convergence rate in space dimension.
	// ####################################

	// Crank-Nicolson. cos(pi * x)
	std::vector<double> space1 = heat_equation_1d(
		"space",
		31,
		"hyperbolic",
		101,
		0.03,
		"exponential",
		"d2dx2::nonuniform::c2b0",
		"cos(pi*x)",
		20,
		0.5);

	// Maximum norm.
	EXPECT_NEAR(space1[0], 2.0, 0.008);

	// L1 function norm.
	EXPECT_NEAR(space1[3], 2.0, 0.011);

	// Crank-Nicolson. cos(3 * pi * x)
	std::vector<double> space2 = heat_equation_1d(
		"space",
		51,
		"hyperbolic",
		201,
		0.03,
		"uniform",
		"d2dx2::nonuniform::c2b0",
		"cos(3*pi*x)",
		20,
		0.5);

	// Maximum norm.
	EXPECT_NEAR(space2[0], 2.0, 0.009);

	// L1 function norm.
	EXPECT_NEAR(space2[3], 2.0, 0.009);

	// Fully implicit.
	std::vector<double> space3 = heat_equation_1d(
		"space",
		31,
		"hyperbolic",
		10001,
		0.03,
		"uniform",
		"d2dx2::nonuniform::c2b0",
		"cos(pi*x)",
		20,
		1.0);

	// Maximum norm.
	EXPECT_NEAR(space3[0], 2.0, 0.028);

	// L1 function norm.
	EXPECT_NEAR(space3[3], 2.0, 0.043);

	// ###################################
	// Convergence rate in time dimension.
	// ###################################

	// Crank-Nicolson. 
	std::vector<double> time1 = heat_equation_1d(
		"time",
		5001,
		"uniform",
		11,
		0.03,
		"exponential",
		"d2dx2::uniform::c2b0",
		"cos(pi*x)",
		20,
		0.5);

	// Maximum norm.
	EXPECT_NEAR(time1[0], 2.0, 0.021);

	// L1 function norm.
	EXPECT_NEAR(time1[3], 2.0, 0.021);

	// Crank-Nicolson. 
	std::vector<double> time2 = heat_equation_1d(
		"time",
		5001,
		"uniform",
		11,
		0.03,
		"uniform",
		"d2dx2::uniform::c2b0",
		"cos(3*pi*x)",
		20,
		0.5);

	// Maximum norm.
	EXPECT_NEAR(time2[0], 2.0, 0.004);

	// L1 function norm.
	EXPECT_NEAR(time2[3], 2.0, 0.004);

	// Fully implicit.
	std::vector<double> time3 = heat_equation_1d(
		"time",
		5001,
		"uniform",
		11,
		0.03,
		"uniform",
		"d2dx2::uniform::c2b0",
		"cos(pi*x)",
		20,
		1.0);

	// Maximum norm.
	EXPECT_NEAR(time3[0], 1.0, 0.057);

	// L1 function norm.
	EXPECT_NEAR(time3[3], 1.0, 0.057);

}


TEST(PentaDiagonalSolver, HeatEquation1D) {

	// Crank-Nicolson. cos(pi * x)
	std::vector<double> space1 = heat_equation_1d(
		"space",
		21,
		"uniform",
		101,
		0.03,
		"uniform",
		"d2dx2::uniform::c4b0",
		"cos(pi*x)",
		11,
		0.5);

	const int n_points = 11;
	std::vector<double> grid_new = grid::uniform(-0.5, 0.5, n_points);
	const double dx = grid_new[1] - grid_new[0];
	PentaDiagonal deriv = d2dx2::uniform::c4b0(n_points, dx);
	print_matrix(deriv);

	PentaDiagonal deriv_non = d2dx2::nonuniform::c4b0(n_points, grid_new);
	print_matrix(deriv_non);

	// Crank-Nicolson. cos(pi * x)
	std::vector<double> space2 = heat_equation_1d(
		"space",
		21,
		"hyperbolic",
		101,
		0.03,
		"uniform",
		"d2dx2::nonuniform::c4b0",
		"cos(pi*x)",
		11,
		0.5);

}
