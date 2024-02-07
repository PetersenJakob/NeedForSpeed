#include "pch.h"


std::vector<double> test_fd_approximation(
	const int function_type,
	const int derivative_order,
	const std::string fd_deriv_type,
	const int n_points_start,
	const int step_size,
	const std::string grid_type,
	bool show_output_basic = false,
	bool show_output_all = false) {
	
	std::vector<double> dx_vec;
	std::vector<double> max_norm;
	std::vector<double> l1_vec_norm;
	std::vector<double> l2_vec_norm;
	std::vector<double> l1_func_norm;
	std::vector<double> l2_func_norm;

	// TODO: Reduce to 50 or 30
	for (int i = 0; i != 100; ++i) {

		// Number of grid points.
		const int n_points = n_points_start + step_size * i;

		// Grid. TODO: The grid is chosen such that the functions are NOT zero at the boundary.
		std::vector<double> grid;

		if (grid_type == "equidistant")
		{
			grid = grid::equidistant(-0.4, 0.4, n_points);
		}
		else if (grid_type == "exponential")
		{
			grid = grid::exponential(-0.4, 0.4, n_points);
		}
		else if (grid_type == "hyperbolic")
		{
			grid = grid::hyperbolic(-0.4, 0.4, n_points);
		}
		else {
			throw std::invalid_argument("grid_type unknown.");
		}

		// (Average) Grid spacing.
		double dx_sum = 0.0;
		for (int i = 0; i != grid.size() - 1; ++i) {
			dx_sum += grid[i + 1] - grid[i];
		}

		const double dx = dx_sum / (grid.size() - 1);
		dx_vec.push_back(dx);

		// Function.
		const std::vector<double> func = test_util::test_function(grid, function_type, 0);

		// Derivative (analytical expression).
		const std::vector<double> deriv = test_util::test_function(grid, function_type, derivative_order);

		// Derivative (finite difference approximaton).
		std::vector<double> deriv_fd;

		if (fd_deriv_type == "d1dx1::equidistant::c2b1") {
			TriDiagonal dndxn = d1dx1::equidistant::c2b1(n_points, dx);
			deriv_fd = dndxn.mat_vec_prod(func);
		}

		else if (fd_deriv_type == "d1dx1::nonequidistant::c2b1") {
			TriDiagonal dndxn = d1dx1::nonequidistant::c2b1(n_points, grid);
			deriv_fd = dndxn.mat_vec_prod(func);
		}

		else if (fd_deriv_type == "d1dx1::equidistant::c2b2") {
			TriDiagonal dndxn = d1dx1::equidistant::c2b2(n_points, dx);
			deriv_fd = dndxn.mat_vec_prod(func);
		}

		else if (fd_deriv_type == "d1dx1::nonequidistant::c2b2") {
			TriDiagonal dndxn = d1dx1::nonequidistant::c2b2(n_points, grid);
			deriv_fd = dndxn.mat_vec_prod(func);
		}

		else if (fd_deriv_type == "d1dx1::equidistant::c4b2") {
			PentaDiagonal dndxn = d1dx1::equidistant::c4b2(n_points, dx);
			deriv_fd = dndxn.mat_vec_prod(func);
		}

		else if (fd_deriv_type == "d1dx1::equidistant::c4b4") {
			PentaDiagonal dndxn = d1dx1::equidistant::c4b4(n_points, dx);
			deriv_fd = dndxn.mat_vec_prod(func);
		}

		else if (fd_deriv_type == "d1dx1::nonequidistant::c4b2") {
			PentaDiagonal dndxn = d1dx1::nonequidistant::c4b2(n_points, grid);
			deriv_fd = dndxn.mat_vec_prod(func);
		}

		else if (fd_deriv_type == "d2dx2::equidistant::c2b1") {
			TriDiagonal dndxn = d2dx2::equidistant::c2b1(n_points, dx);
			deriv_fd = dndxn.mat_vec_prod(func);
		}
		else if (fd_deriv_type == "d2dx2::equidistant::c2b2") {
			TriDiagonal dndxn = d2dx2::equidistant::c2b2(n_points, dx);
			deriv_fd = dndxn.mat_vec_prod(func);
		}
		else if (fd_deriv_type == "d2dx2::equidistant::c4b4") {
			PentaDiagonal dndxn = d2dx2::equidistant::c4b4(n_points, dx);
			deriv_fd = dndxn.mat_vec_prod(func);
		}
		else {
			throw std::invalid_argument("FD derivative unknown.");
		}

		if (show_output_all) {
			test_util::print_test_results(grid, func, deriv, deriv_fd);
		}

		// Difference vector.
		std::vector<double> diff = test_util::vector_diff(deriv, deriv_fd);

		max_norm.push_back(test_util::max_norm(diff));

		l1_vec_norm.push_back(test_util::l1_vector_norm(diff));

		l2_vec_norm.push_back(test_util::l2_vector_norm(diff));

		l1_func_norm.push_back(test_util::l1_function_norm(dx, diff));

		l2_func_norm.push_back(test_util::l2_function_norm(dx, diff));

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

		std::ofstream myfile("slr.txt");
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

	return result;

}


TEST(FirstOrderDerivative, EXPc2b1) {

	std::vector<double> slope = test_fd_approximation(0, 1, "d1dx1::equidistant::c2b1", 51, 5, "equidistant", false, false);

	// Maximum norm.
	EXPECT_NEAR(slope[0], 1.0, 0.004);

	// L1 function norm.
	EXPECT_NEAR(slope[3], 2.0, 0.008);

}


TEST(FirstOrderDerivative, COSc2b1) {

	std::vector<double> slope = test_fd_approximation(1, 1, "d1dx1::equidistant::c2b1", 51, 5, "equidistant", true, false);

	// Maximum norm.
	EXPECT_NEAR(slope[0], 1.0, 0.015);

	// L1 function norm.
	EXPECT_NEAR(slope[3], 2.0, 0.001);

}


TEST(FirstOrderDerivative, SUMc2b1) {

	std::vector<double> slope = test_fd_approximation(2, 1, "d1dx1::equidistant::c2b1", 51, 5, "equidistant", true, false);

	// Maximum norm.
	EXPECT_NEAR(slope[0], 1.0, 0.013);

	// L1 function norm.
	EXPECT_NEAR(slope[3], 2.0, 0.017);

}


TEST(FirstOrderDerivative, EXPc2b2) {


	std::vector<double> slope = test_fd_approximation(0, 1, "d1dx1::equidistant::c2b2", 51, 5, "equidistant", true, false);

	// Maximum norm.
	EXPECT_NEAR(slope[0], 2.0, 0.008);

	// L1 function norm.
	EXPECT_NEAR(slope[3], 2.0, 0.008);

}


TEST(FirstOrderDerivative, COSc2b2) {

	std::vector<double> slope = test_fd_approximation(1, 1, "d1dx1::equidistant::c2b2", 51, 5, "equidistant", true, false);

	// Maximum norm.
	EXPECT_NEAR(slope[0], 2.0, 0.004);

	// L1 function norm.
	EXPECT_NEAR(slope[3], 2.0, 0.011);

}


TEST(FirstOrderDerivative, SUMc2b2) {

	std::vector<double> slope = test_fd_approximation(2, 1, "d1dx1::equidistant::c2b2", 51, 5, "equidistant", true, false);

	// Maximum norm.
	EXPECT_NEAR(slope[0], 2.0, 0.006);

	// L1 function norm.
	EXPECT_NEAR(slope[3], 2.0, 0.011);

}


TEST(FirstOrderDerivative, EXPc4b4) {

	std::vector<double> slope = test_fd_approximation(0, 1, "d1dx1::equidistant::c4b4", 51, 5, "equidistant", true, false);

	// Maximum norm.
	EXPECT_NEAR(slope[0], 4.0, 0.016);

	// L1 function norm.
	EXPECT_NEAR(slope[3], 4.0, 0.090);

}


TEST(FirstOrderDerivative, COSc4b4) {

	std::vector<double> slope = test_fd_approximation(1, 1, "d1dx1::equidistant::c4b4", 51, 5, "equidistant", true, false);

	// Maximum norm.
	EXPECT_NEAR(slope[0], 4.0, 0.010);

	// L1 function norm.
	EXPECT_NEAR(slope[3], 4.0, 0.124);

}


TEST(FirstOrderDerivative, SUMc4b4) {

	std::vector<double> slope = test_fd_approximation(2, 1, "d1dx1::equidistant::c4b4", 51, 5, "equidistant", true, false);

	// Maximum norm.
	EXPECT_NEAR(slope[0], 4.0, 0.009);

	// L1 function norm.
	EXPECT_NEAR(slope[3], 4.0, 0.121);

}


TEST(SecondOrderDerivative, EXPc2b1) {

	std::vector<double> slope = test_fd_approximation(0, 2, "d2dx2::equidistant::c2b1", 51, 5, "equidistant", true, false);

	// Maximum norm.
	EXPECT_NEAR(slope[0], 1.0, 0.006);

	// L1 function norm.
	EXPECT_NEAR(slope[3], 2.0, 0.006);

}


TEST(SecondOrderDerivative, COSc2b1) {

	std::vector<double> slope = test_fd_approximation(1, 2, "d2dx2::equidistant::c2b1", 51, 5, "equidistant", true, false);

	// Maximum norm.
	EXPECT_NEAR(slope[0], 1.0, 0.003);

	// L1 function norm.
	EXPECT_NEAR(slope[3], 2.0, 0.004);

}


TEST(SecondOrderDerivative, SUMc2b1) {

	std::vector<double> slope = test_fd_approximation(2, 2, "d2dx2::equidistant::c2b1", 51, 5, "equidistant", true, false);

	// Maximum norm.
	EXPECT_NEAR(slope[0], 1.0, 0.004);

	// L1 function norm.
	EXPECT_NEAR(slope[3], 2.0, 0.005);

}


TEST(SecondOrderDerivative, EXPc2b2) {

	std::vector<double> slope = test_fd_approximation(0, 2, "d2dx2::equidistant::c2b2", 51, 5, "equidistant", true, false);

	// Maximum norm.
	EXPECT_NEAR(slope[0], 2.0, 0.011);

	// L1 function norm.
	EXPECT_NEAR(slope[3], 2.0, 0.052);

}


TEST(SecondOrderDerivative, COSc2b2) {

	std::vector<double> slope = test_fd_approximation(1, 2, "d2dx2::equidistant::c2b2", 51, 5, "equidistant", true, false);

	// Maximum norm.
	EXPECT_NEAR(slope[0], 2.0, 0.047);

	// L1 function norm.
	EXPECT_NEAR(slope[3], 2.0, 0.022);

}


TEST(SecondOrderDerivative, SUMc2b2) {

	std::vector<double> slope = test_fd_approximation(2, 2, "d2dx2::equidistant::c2b2", 51, 5, "equidistant", true, false);

	// Maximum norm.
	EXPECT_NEAR(slope[0], 2.0, 0.017);

	// L1 function norm.
	EXPECT_NEAR(slope[3], 2.0, 0.028);

}


TEST(SecondOrderDerivative, EXPc4b4) {

	std::vector<double> slope = test_fd_approximation(0, 2, "d2dx2::equidistant::c4b4", 51, 2, "equidistant", true, false);

	// Maximum norm.
	EXPECT_NEAR(slope[0], 4.0, 0.063);

	// L1 function norm.
	EXPECT_NEAR(slope[3], 4.0, 0.600);

}


TEST(SecondOrderDerivative, COSc4b4) {

	std::vector<double> slope = test_fd_approximation(1, 2, "d2dx2::equidistant::c4b4", 51, 2, "equidistant", true, false);

	// Maximum norm.
	EXPECT_NEAR(slope[0], 4.0, 0.166);

	// L1 function norm.
	EXPECT_NEAR(slope[3], 4.0, 0.528);

}


TEST(SecondOrderDerivative, SUMc4b4) {

	std::vector<double> slope = test_fd_approximation(2, 2, "d2dx2::equidistant::c4b4", 51, 2, "equidistant", true, false);

	// Maximum norm.
	EXPECT_NEAR(slope[0], 4.0, 0.181);

	// L1 function norm.
	EXPECT_NEAR(slope[3], 4.0, 0.516);

}



// TODO: Combine equidistant and non-equidistant tests...
TEST(FirstOrderDerivativeNonequidistant, EXPc2b1) {

	// Compare FD representations based on equidistant and non-equidistant grids.
	
	const int n_points = 21;

	// Grid.
	const std::vector<double> grid_eq = grid::equidistant(-0.4, 0.4, n_points);

	const std::vector<double> grid_exp = grid::exponential(-0.4, 0.4, n_points);

	const std::vector<double> grid_hyper = grid::hyperbolic(-0.4, 0.4, n_points);

#if false
	std::ofstream myfile("grid.csv");
	myfile << std::scientific << std::setprecision(12);
	for (int i = 0; i != grid_eq.size(); ++i) {

		myfile 
			<< std::setw(22) << grid_eq[i] << ","
			<< std::setw(22) << grid_exp[i] << ","
			<< std::setw(22) << grid_hyper[i] << std::endl

	}
	myfile << std::endl;
	myfile.close();
#endif

	// Grid spacing.
	const double dx = grid_eq[1] - grid_eq[0];

	TriDiagonal d1dx1_eq = d1dx1::equidistant::c2b1(n_points, dx);
	TriDiagonal d1dx1_neq = d1dx1::nonequidistant::c2b1(n_points, grid_eq);

//	print_matrix(d1dx1_eq);
//	print_matrix(d1dx1_neq);

	EXPECT_TRUE(d1dx1_eq == d1dx1_neq);

	std::vector<double> slope_eq = test_fd_approximation(0, 1, "d1dx1::equidistant::c2b1", 51, 5, "equidistant", true, false);

	std::vector<double> slope_exp = test_fd_approximation(0, 1, "d1dx1::nonequidistant::c2b1", 51, 5, "exponential", true, false);

	std::vector<double> slope_hyper = test_fd_approximation(0, 1, "d1dx1::nonequidistant::c2b1", 51, 5, "hyperbolic", true, false);

	// Maximum norm.
	EXPECT_NEAR(slope_eq[0], 1.0, 0.004);

	EXPECT_NEAR(slope_exp[0], 1.0, 0.013);

	EXPECT_NEAR(slope_hyper[0], 1.0, 0.021);

	// L1 function norm.
	EXPECT_NEAR(slope_eq[3], 2.0, 0.008);

	EXPECT_NEAR(slope_exp[3], 2.0, 0.025);

	EXPECT_NEAR(slope_hyper[3], 2.0, 0.029);

}


TEST(FirstOrderDerivativeNonequidistant, EXPc2b2) {

	// Compare FD representations based on equidistant and non-equidistant grids.

	const int n_points = 21;

	// Grid.
	const std::vector<double> grid_eq = grid::equidistant(-0.4, 0.4, n_points);

	const std::vector<double> grid_exp = grid::exponential(-0.4, 0.4, n_points);

	const std::vector<double> grid_hyper = grid::hyperbolic(-0.4, 0.4, n_points);

	std::ofstream myfile("grid.csv");
	myfile << std::scientific << std::setprecision(12);
	for (int i = 0; i != grid_eq.size(); ++i) {

		myfile
			<< std::setw(22) << grid_eq[i] << ","
			<< std::setw(22) << grid_exp[i] << ","
			<< std::setw(22) << grid_hyper[i] << std::endl;

	}
	myfile << std::endl;
	myfile.close();

	// Grid spacing.
	const double dx = grid_eq[1] - grid_eq[0];

	TriDiagonal d1dx1_eq = d1dx1::equidistant::c2b2(n_points, dx);
	TriDiagonal d1dx1_neq = d1dx1::nonequidistant::c2b2(n_points, grid_eq);

	//	print_matrix(d1dx1_eq);
	//	print_matrix(d1dx1_neq);

	EXPECT_TRUE(d1dx1_eq == d1dx1_neq);

	std::vector<double> slope_eq = test_fd_approximation(0, 1, "d1dx1::equidistant::c2b2", 51, 5, "equidistant", true, false);

	std::vector<double> slope_exp = test_fd_approximation(0, 1, "d1dx1::nonequidistant::c2b2", 51, 5, "exponential", true, false);

	std::vector<double> slope_hyper = test_fd_approximation(0, 1, "d1dx1::nonequidistant::c2b2", 51, 5, "hyperbolic", true, false);

	// Maximum norm.
	EXPECT_NEAR(slope_eq[0], 2.0, 0.008);

	EXPECT_NEAR(slope_exp[0], 2.0, 0.034);

	EXPECT_NEAR(slope_hyper[0], 2.0, 0.056);

	// L1 function norm.
	EXPECT_NEAR(slope_eq[3], 2.0, 0.008);

	EXPECT_NEAR(slope_exp[3], 2.0, 0.019);

	EXPECT_NEAR(slope_hyper[3], 2.0, 0.027);

}


TEST(FirstOrderDerivativeNonequidistant, EXPc4b2) {

	// Compare FD representations based on equidistant and non-equidistant grids.

	const int n_points = 21;

	// Grid.
	const std::vector<double> grid_eq = grid::equidistant(-0.4, 0.4, n_points);

	const std::vector<double> grid_exp = grid::exponential(-0.4, 0.4, n_points);

	const std::vector<double> grid_hyper = grid::hyperbolic(-0.4, 0.4, n_points);


	// Grid spacing.
	const double dx = grid_eq[1] - grid_eq[0];

	PentaDiagonal d1dx1_eq = d1dx1::equidistant::c4b2(n_points, dx);
	PentaDiagonal d1dx1_neq = d1dx1::nonequidistant::c4b2(n_points, grid_eq);

//	print_matrix(d1dx1_eq);
//	print_matrix(d1dx1_neq);

	EXPECT_TRUE(d1dx1_eq == d1dx1_neq);

	std::vector<double> slope_eq = test_fd_approximation(0, 1, "d1dx1::equidistant::c4b2", 51, 5, "equidistant", true, false);

	std::vector<double> slope_exp = test_fd_approximation(0, 1, "d1dx1::nonequidistant::c4b2", 51, 5, "exponential", true, false);

	std::vector<double> slope_hyper = test_fd_approximation(0, 1, "d1dx1::nonequidistant::c4b2", 51, 5, "hyperbolic", true, false);

	// Maximum norm.
	EXPECT_NEAR(slope_eq[0], 2.0, 0.008);

	EXPECT_NEAR(slope_exp[0], 2.0, 0.034);

	EXPECT_NEAR(slope_hyper[0], 2.0, 0.056);

	// L1 function norm.
	EXPECT_NEAR(slope_eq[3], 3.0, 0.009);

	EXPECT_NEAR(slope_exp[3], 3.0, 0.059);

	EXPECT_NEAR(slope_hyper[3], 3.0, 0.084);

}
