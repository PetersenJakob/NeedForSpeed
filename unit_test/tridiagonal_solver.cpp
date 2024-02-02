#include "pch.h"


TEST(TriDiagonalSolver, HeatEquation1D) {

	bool show_output_basic = false;
	bool show_output_all = false;

	std::vector<double> dx_vec;
	std::vector<double> max_norm;
	std::vector<double> l2_norm;

	// Time interval.
	const double time_interval = 0.1; // 0.1;

	// Number of time steps.
	int n_steps = 101;

	// Number of grid points.
	int n_points_start = 21;
	int n_points = 0;

	for (int i = 0; i != 15; ++i) {

		// Step size in time.
		double dt = time_interval / (n_steps - 1);

		// Number of grid points.
		n_points = n_points_start + 10 * i;

		// Grid.
		const std::vector<double> grid = grid_equidistant(-0.5, 0.5, n_points);

		// Grid spacing.
		const double dx = grid[1] - grid[0];
		dx_vec.push_back(dx);

		// Initial condition.
		std::vector<double> func(n_points, 0.0);
		for (int i = 0; i != n_points; ++i) {
			func[i] = cos(M_PI * grid[i]);
		}

		// Second order differential operator.
		TriDiagonal deriv_operator = d2dx2::c2b1(n_points, dx);


		std::vector<double> coefficients(3, 0.0);
		boundary<TriDiagonal>(0, pow(dx, 2.0), coefficients, deriv_operator);
		boundary<TriDiagonal>(1, pow(dx, 2.0), coefficients, deriv_operator);

//		print_matrix(deriv_operator);

#if false
		PentaDiagonal deriv_operator_p = d2dx2::c4b4(n_points, dx);

		std::vector<double> coefficients_p(6, 0.0);
		boundary<PentaDiagonal>(0, pow(dx, 2.0), coef2::f1, deriv_operator_p);
		boundary<PentaDiagonal>(1, pow(dx, 2.0), coefficients_p, deriv_operator_p);
		boundary<PentaDiagonal>(2, pow(dx, 2.0), coefficients_p, deriv_operator_p);
		boundary<PentaDiagonal>(3, pow(dx, 2.0), coefficients_p, deriv_operator_p);

//		print_matrix(deriv_operator_p); 
#endif

		// Theta parameter.
		const double theta = 0.5;

		// LHS operator.
		TriDiagonal lhs = deriv_operator;
		lhs.scalar_prod(-theta * dt);
		lhs.add_diagonal(1.0);

		// RHS operator.
		TriDiagonal rhs = deriv_operator;
		rhs.scalar_prod((1.0 - theta) * dt);
		rhs.add_diagonal(1.0);


#if false
		PentaDiagonal lhs_p = deriv_operator_p;
		lhs_p.scalar_prod(-theta * dt);
		lhs_p.add_diagonal(1.0);

		PentaDiagonal rhs_p = deriv_operator_p;
		rhs_p.scalar_prod((1.0 - theta) * dt);
		rhs_p.add_diagonal(1.0);

//		print_matrix(lhs_p);
//		print_matrix(rhs_p);
#endif

		// FD solution.
		std::vector<double> solution_fd = func;

		for (int i = 0; i != n_steps; ++i) {


			solution_fd = rhs.mat_vec_prod(solution_fd);
			tri_solver(lhs, solution_fd);

#if false
			solution_fd = rhs_p.mat_vec_prod(solution_fd);
			penta_solver(lhs_p, solution_fd);
#endif

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

		std::vector<double> diff = test_util::vector_diff(analytical, solution_fd);
		max_norm.push_back(test_util::max_norm(diff));
		l2_norm.push_back(test_util::l2_norm(dx, diff));

	}

	std::vector<double> tmp(3, 0.0);
	std::vector<std::vector<double>> ratio(2, tmp);

	ratio[0][0] = max_norm[0] / max_norm[2];
	ratio[0][1] = max_norm[2] / max_norm[6];
	ratio[0][2] = max_norm[6] / max_norm[14];

	ratio[1][0] = l2_norm[0] / l2_norm[2];
	ratio[1][1] = l2_norm[2] / l2_norm[6];
	ratio[1][2] = l2_norm[6] / l2_norm[14];

	if (show_output_basic) {
		std::cout << std::scientific << std::setprecision(5);
		std::cout << "    dx     max-norm       l2-norm" << std::endl;
		for (int i = 0; i != max_norm.size(); ++i) {
			std::cout
				<< std::setw(14) << dx_vec[i]
				<< std::setw(14) << max_norm[i]
				<< std::setw(14) << l2_norm[i]
				<< std::endl;
		}
		std::cout << std::endl;


		std::cout
			<< "Max-norm:" << std::endl
			<< ratio[0][0] << std::endl
			<< ratio[0][1] << std::endl
			<< ratio[0][2] << std::endl << std::endl
			<< "l2-norm:" << std::endl
			<< ratio[1][0] << std::endl
			<< ratio[1][1] << std::endl
			<< ratio[1][2] << std::endl << std::endl;
	}

	for (int i = 0; i != ratio[0].size(); ++i) {
		EXPECT_NEAR(ratio[0][i], 4.5, 0.6);
	}

	for (int i = 0; i != ratio[1].size(); ++i) {
		EXPECT_NEAR(ratio[1][i], 4.5, 0.6);
	}

}


TEST(TriDiagonalSolver, BlackScholes1D) {

	bool show_output_basic = false;
	bool show_output_all = false;

	std::vector<double> dx_vec;
	std::vector<double> max_norm;
	std::vector<double> l2_norm;

	// Time interval.
	const double time_interval = 0.1;

	// Number of time steps.
	int n_steps = 101;

	// Number of grid points.
	int n_points = 21;

	for (int i = 0; i != 15; ++i) {

		// Step size in time.
		double dt = time_interval / (n_steps - 1);

		// Number of grid points.
		n_points = 21 + 10 * i;

		// Grid.
		const std::vector<double> grid = grid_equidistant(-0.5, 0.5, n_points);

		// Grid spacing.
		const double dx = grid[1] - grid[0];
		dx_vec.push_back(dx);

		// Initial condition.
		std::vector<double> func(n_points, 0.0);
		for (int i = 0; i != n_points; ++i) {
			func[i] = cos(M_PI * grid[i]);
		}

		// Second order differential operator.
		TriDiagonal deriv_operator = d2dx2::c2b1(n_points, dx);


		std::vector<double> coefficients(3, 0.0);
		boundary<TriDiagonal>(0, pow(dx, 2.0), coefficients, deriv_operator);
		boundary<TriDiagonal>(1, pow(dx, 2.0), coefficients, deriv_operator);


		// Theta parameter.
		const double theta = 0.5;

		// LHS operator.
		TriDiagonal lhs = deriv_operator;
		lhs.scalar_prod(-theta * dt);
		lhs.add_diagonal(1.0);

		// RHS operator.
		TriDiagonal rhs = deriv_operator;
		rhs.scalar_prod((1.0 - theta) * dt);
		rhs.add_diagonal(1.0);

		// FD solution.
		std::vector<double> solution_fd = func;

		for (int i = 0; i != n_steps; ++i) {

			solution_fd = rhs.mat_vec_prod(solution_fd);
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

		std::vector<double> diff = test_util::vector_diff(analytical, solution_fd);
		max_norm.push_back(test_util::max_norm(diff));
		l2_norm.push_back(test_util::l2_norm(dx, diff));

	}

	std::vector<double> tmp(3, 0.0);
	std::vector<std::vector<double>> ratio(2, tmp);

	ratio[0][0] = max_norm[0] / max_norm[2];
	ratio[0][1] = max_norm[2] / max_norm[6];
	ratio[0][2] = max_norm[6] / max_norm[14];

	ratio[1][0] = l2_norm[0] / l2_norm[2];
	ratio[1][1] = l2_norm[2] / l2_norm[6];
	ratio[1][2] = l2_norm[6] / l2_norm[14];

	if (show_output_basic) {
		std::cout << std::scientific << std::setprecision(5);
		std::cout << "    dx     max-norm       l2-norm" << std::endl;
		for (int i = 0; i != max_norm.size(); ++i) {
			std::cout
				<< std::setw(14) << dx_vec[i]
				<< std::setw(14) << max_norm[i]
				<< std::setw(14) << l2_norm[i]
				<< std::endl;
		}
		std::cout << std::endl;


		std::cout
			<< "Max-norm:" << std::endl
			<< ratio[0][0] << std::endl
			<< ratio[0][1] << std::endl
			<< ratio[0][2] << std::endl << std::endl
			<< "l2-norm:" << std::endl
			<< ratio[1][0] << std::endl
			<< ratio[1][1] << std::endl
			<< ratio[1][2] << std::endl << std::endl;
	}

	for (int i = 0; i != ratio[0].size(); ++i) {
		EXPECT_NEAR(ratio[0][i], 4.5, 0.6);
	}

	for (int i = 0; i != ratio[1].size(); ++i) {
		EXPECT_NEAR(ratio[1][i], 4.5, 0.6);
	}

}
