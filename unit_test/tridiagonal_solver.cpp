#include "pch.h"


TEST(TriDiagonalSolver, HeatEquation1D) {

	bool show_output_basic = true;
	bool show_output_all = false;

	std::vector<double> dx_vec;
	std::vector<double> max_norm;
	std::vector<double> l1_vec_norm;
	std::vector<double> l2_vec_norm;
	std::vector<double> l1_func_norm;
	std::vector<double> l2_func_norm;

	// Time interval.
	const double time_interval = 0.04;

	// Number of time steps.
	int n_steps = 251;

	// Number of grid points.
	int n_points_start = 51;
	int n_points = 0;

	// Theta parameter.
	const double theta = 0.5;

	for (int i = 0; i != 50; ++i) {

		// Step size in time.
		double dt = time_interval / (n_steps - 1);

		// Number of grid points.
		n_points = n_points_start + 5 * i;

		// Grid. TODO: Chosen such that d2dx2 = 0 for analytical solution at boundary!
//		const std::vector<double> grid = grid::equidistant(-0.5, 0.5, n_points);
		const std::vector<double> grid = grid::hyperbolic(-0.5, 0.5, n_points);

		// (Average) Grid spacing.
		double dx_sum = 0.0;
		for (int i = 0; i != grid.size() - 1; ++i) {
			dx_sum += grid[i + 1] - grid[i];
		}
		const double dx = dx_sum / (grid.size() - 1);
		dx_vec.push_back(dx);


		// Initial condition.
		std::vector<double> func(n_points, 0.0);
		for (int i = 0; i != n_points; ++i) {
			func[i] = cos(M_PI * grid[i]);
		}


		// Second order differential operator.
//		TriDiagonal deriv_operator = d2dx2::equidistant::c2b0(n_points, dx);
		TriDiagonal deriv_operator = d2dx2::nonequidistant::c2b0(n_points, grid);

		// Identity operator.
		TriDiagonal iden = identity::tri(n_points, 2);

		// LHS operator.
		TriDiagonal lhs = deriv_operator;
		lhs *= -theta * dt;
		lhs += iden;

		// RHS operator.
		TriDiagonal rhs = deriv_operator;
		rhs *= (1.0 - theta) * dt;
		rhs += iden;

		// FD solution.
		std::vector<double> solution_fd = func;

		for (int i = 0; i != n_steps; ++i) {

			solution_fd = rhs * solution_fd;

			tri_solver(lhs, solution_fd);

		}

		// Analytical solution.
		std::vector<double> analytical(n_points, 0.0);
		for (int i = 0; i != n_points; ++i) {
			analytical[i] = exp(-pow(M_PI, 2.0) * n_steps * dt) * func[i];
		}

		if (show_output_all) {
			std::cout << std::scientific << std::setprecision(5);
			for (int i = 0; i != n_points; ++i) {
				std::cout
					<< std::setw(3) << i
					<< std::setw(14) << grid[i]
					<< std::setw(14) << func[i]
					<< std::setw(14) << analytical[i]
					<< std::setw(14) << solution_fd[i]
					<< std::setw(14) << abs(analytical[i] - solution_fd[i])
					<< std::endl;

			}
			std::cout << std::endl;
		}

		// Difference vector.
		std::vector<double> diff = test_util::vector_diff(analytical, solution_fd);

		max_norm.push_back(test_util::max_norm(diff));

		l1_vec_norm.push_back(test_util::l1_vector_norm(diff));

		l2_vec_norm.push_back(test_util::l2_vector_norm(diff));

		l1_func_norm.push_back(test_util::l1_function_norm(grid, diff));

		l2_func_norm.push_back(test_util::l2_function_norm(grid, diff));

	}

	std::vector<double> dx_vec_log(dx_vec.size(), 0.0);
	std::vector<double> max_norm_log(max_norm.size(), 0.0);
	std::vector<double> l1_vec_norm_log(l1_vec_norm.size(), 0.0);
	std::vector<double> l2_vec_norm_log(l2_vec_norm.size(), 0.0);
	std::vector<double> l1_func_norm_log(l1_func_norm.size(), 0.0);
	std::vector<double> l2_func_norm_log(l2_func_norm.size(), 0.0);

	for (int i = 0; i != dx_vec.size(); ++i) {

		dx_vec_log[i] = log(dx_vec[i]);
		max_norm_log[i] = log(max_norm[i]);
		l1_vec_norm_log[i] = log(l1_vec_norm[i]);
		l2_vec_norm_log[i] = log(l2_vec_norm[i]);
		l1_func_norm_log[i] = log(l1_func_norm[i]);
		l2_func_norm_log[i] = log(l2_func_norm[i]);

	}

	std::vector<double> slr_max = test_util::slr(dx_vec_log, max_norm_log);
	std::vector<double> slr_l1_vec = test_util::slr(dx_vec_log, l1_vec_norm_log);
	std::vector<double> slr_l2_vec = test_util::slr(dx_vec_log, l2_vec_norm_log);
	std::vector<double> slr_l1_func = test_util::slr(dx_vec_log, l1_func_norm_log);
	std::vector<double> slr_l2_func = test_util::slr(dx_vec_log, l2_func_norm_log);

	std::vector<double> result{ slr_max[0],
		slr_l1_vec[0], slr_l2_vec[0],
		slr_l1_func[0], slr_l2_func[0] };

	if (show_output_basic) {

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

		std::ofstream myfile("slr.csv");
		myfile << std::scientific << std::setprecision(12);
		for (int i = 0; i != dx_vec.size(); ++i) {

			myfile
				<< std::setw(22) << dx_vec[i] << ","
				<< std::setw(22) << dx_vec_log[i] << ","
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
		myfile.close();

	}

	// Maximum norm.
	EXPECT_NEAR(result[0], 2.0, 0.007);

	// L1 function norm.
	EXPECT_NEAR(result[3], 2.0, 0.007);

}
