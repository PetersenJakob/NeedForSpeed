// Unit tests of finite difference approximations!!!!

#include "pch.h"


std::vector<double> test_fd_approximation(
	const int function_type,
	const int derivative_order,
	const std::string fd_deriv_type,
	bool show_output_basic = false,
	bool show_output_all = false) {
	
	std::vector<double> dx_vec;
	std::vector<double> max_norm;
	std::vector<double> l2_norm;

	for (int i = 0; i != 50; ++i) {

		// Number of grid points.
		const int n_points = 101 + 2 * i;

		// Grid.
		const std::vector<double> grid = grid_equidistant(-0.2, 0.2, n_points);

		// Grid spacing.
		const double dx = grid[1] - grid[0];
		dx_vec.push_back(dx);

		// Function.
		const std::vector<double> func = test_util::test_function(grid, function_type, 0);

		// Derivative (analytical expression).
		const std::vector<double> deriv = test_util::test_function(grid, function_type, derivative_order);

		// Derivative (finite difference approximaton).

		std::vector<double> deriv_fd;

		if (fd_deriv_type == "d1dx1::c2b1") {
			TriDiagonal dndxn = d1dx1::equidistant::c2b1(n_points, dx);
			deriv_fd = dndxn.mat_vec_prod(func);
		}
		else if (fd_deriv_type == "d1dx1::c2b2") {
			TriDiagonal dndxn = d1dx1::equidistant::c2b2(n_points, dx);

//			print_matrix(dndxn);

			deriv_fd = dndxn.mat_vec_prod(func);
		}
		else if (fd_deriv_type == "d1dx1::c4b4") {
			PentaDiagonal dndxn = d1dx1::equidistant::c4b4(n_points, dx);
			deriv_fd = dndxn.mat_vec_prod(func);
		}
		else if (fd_deriv_type == "d2dx2::c2b1") {
			TriDiagonal dndxn = d2dx2::equidistant::c2b1(n_points, dx);
			deriv_fd = dndxn.mat_vec_prod(func);
		}
		else if (fd_deriv_type == "d2dx2::c2b2") {
			TriDiagonal dndxn = d2dx2::equidistant::c2b2(n_points, dx);
			deriv_fd = dndxn.mat_vec_prod(func);
		}
		else if (fd_deriv_type == "d2dx2::c4b4") {
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




	std::vector<double> dx_vec_log(dx_vec.size(), 0.0);
	std::vector<double> max_norm_log(max_norm.size(), 0.0);
	std::vector<double> l2_norm_log(l2_norm.size(), 0.0);

	for (int i = 0; i != dx_vec.size(); ++i) {

		dx_vec_log[i] = log(dx_vec[i]);
		max_norm_log[i] = log(max_norm[i]);
		l2_norm_log[i] = log(l2_norm[i]);

	}

	std::vector<double> slr_max = test_util::slr(dx_vec_log, max_norm_log);
	std::vector<double> slr_l2 = test_util::slr(dx_vec_log, l2_norm_log);

	std::vector<double> result{ slr_max[0], slr_l2[0] };

	if (show_output_basic) {

		std::cout << std::scientific << std::setprecision(5);

		std::cout 
			<< "SLR max-norm:" << std::endl
			<< std::setw(14) << slr_max[0] 
			<< std::setw(14) << slr_max[1] << std::endl << std::endl;

		std::cout << 
			"SLR l2-norm:" << std::endl
			<< std::setw(14) << slr_l2[0] 
			<< std::setw(14) << slr_l2[1] << std::endl << std::endl;

		std::ofstream myfile("slr.txt");
		myfile << std::scientific << std::setprecision(12);
		for (int i = 0; i != dx_vec.size(); ++i) {

			myfile
				<< std::setw(22) << dx_vec[i] << ","
				<< std::setw(22) << dx_vec_log[i] << ","
				<< std::setw(22) << max_norm[i] << ","
				<< std::setw(22) << max_norm_log[i] << ","
				<< std::setw(22) << l2_norm[i] << ","
				<< std::setw(22) << l2_norm_log[i] << std::endl;

		}
		myfile.close();

	}

	return result;

}


TEST(FirstOrderDerivative, EXPc2b1) {

	std::vector<double> slope = test_fd_approximation(0, 1, "d1dx1::c2b1", true, false);

	// Max norm.
	EXPECT_NEAR(slope[0], 1.0, 0.002);

	// l2-norm.
	EXPECT_NEAR(slope[1], 2.0, 0.003);

}


TEST(FirstOrderDerivative, COSc2b1) {

	std::vector<double> slope = test_fd_approximation(1, 1, "d1dx1::c2b1", true, false);

	// Max norm.
	EXPECT_NEAR(slope[0], 1.0, 0.002);

	// l2-norm.
	EXPECT_NEAR(slope[1], 2.0, 0.001);

}


TEST(FirstOrderDerivative, SUMc2b1) {

	std::vector<double> slope = test_fd_approximation(2, 1, "d1dx1::c2b1", true, false);

	// Max norm.
	EXPECT_NEAR(slope[0], 1.0, 0.002);

	// l2-norm.
	EXPECT_NEAR(slope[1], 2.0, 0.001);

}


TEST(FirstOrderDerivative, EXPc2b2) {


	std::vector<double> slope = test_fd_approximation(0, 1, "d1dx1::c2b2", true, false);

	// Max norm.
	EXPECT_NEAR(slope[0], 2.0, 0.003);

	// l2-norm.
	EXPECT_NEAR(slope[1], 2.8, 0.005);

}


TEST(FirstOrderDerivative, COSc2b2) {

	std::vector<double> slope = test_fd_approximation(1, 1, "d1dx1::c2b2", true, false);

	// Max norm.
	EXPECT_NEAR(slope[0], 2.0, 0.006);

	// l2-norm.
	EXPECT_NEAR(slope[1], 2.0, 0.009);

}


TEST(FirstOrderDerivative, SUMc2b2) {

	std::vector<double> slope = test_fd_approximation(2, 1, "d1dx1::c2b2", true, false);

	// Max norm.
	EXPECT_NEAR(slope[0], 2.0, 0.005);

	// l2-norm.
	EXPECT_NEAR(slope[1], 2.8, 0.008);

}


TEST(FirstOrderDerivative, EXPc4b4) {

	std::vector<double> slope = test_fd_approximation(0, 1, "d1dx1::c4b4", true, false);

	// Max norm.
	EXPECT_NEAR(slope[0], 4.0, 0.01);

	// l2-norm.
	EXPECT_NEAR(slope[1], 4.0, 0.07);

}


TEST(FirstOrderDerivative, COSc4b4) {

	std::vector<double> slope = test_fd_approximation(1, 1, "d1dx1::c4b4", true, false);

	// Max norm.
	EXPECT_NEAR(slope[0], 2.0, 0.02);

	// l2-norm.
	EXPECT_NEAR(slope[1], 2.8, 0.03);

}


TEST(FirstOrderDerivative, SUMc4b4) {

	std::vector<double> slope = test_fd_approximation(2, 1, "d1dx1::c4b4", true, false);

	// Max norm.
	EXPECT_NEAR(slope[0], 2.0, 0.02);

	// l2-norm.
	EXPECT_NEAR(slope[1], 2.8, 0.03);

}


TEST(SecondOrderDerivative, EXPc2b1) {

	std::vector<double> slope = test_fd_approximation(0, 2, "d2dx2::c2b1", true, false);

	// Max norm.
	EXPECT_NEAR(slope[0], 2.0, 0.02);

	// l2-norm.
	EXPECT_NEAR(slope[1], 2.8, 0.03);

}


TEST(SecondOrderDerivative, COSc2b1) {

	std::vector<double> slope = test_fd_approximation(1, 2, "d2dx2::c2b1", true, false);

	// Max norm.
	EXPECT_NEAR(slope[0], 2.0, 0.02);

	// l2-norm.
	EXPECT_NEAR(slope[1], 2.8, 0.03);

}


TEST(SecondOrderDerivative, SUMc2b1) {

	std::vector<double> slope = test_fd_approximation(2, 2, "d2dx2::c2b1", true, false);

	// Max norm.
	EXPECT_NEAR(slope[0], 2.0, 0.02);

	// l2-norm.
	EXPECT_NEAR(slope[1], 2.8, 0.03);

}


TEST(SecondOrderDerivative, EXPc2b2) {

	std::vector<double> slope = test_fd_approximation(0, 2, "d2dx2::c2b2", true, false);

	// Max norm.
	EXPECT_NEAR(slope[0], 2.0, 0.02);

	// l2-norm.
	EXPECT_NEAR(slope[1], 2.8, 0.03);

}


TEST(SecondOrderDerivative, COSc2b2) {

	std::vector<double> slope = test_fd_approximation(1, 2, "d2dx2::c2b2", true, false);

	// Max norm.
	EXPECT_NEAR(slope[0], 2.0, 0.02);

	// l2-norm.
	EXPECT_NEAR(slope[1], 2.8, 0.03);

}


TEST(SecondOrderDerivative, SUMc2b2) {

	std::vector<double> slope = test_fd_approximation(2, 2, "d2dx2::c2b2", true, false);

	// Max norm.
	EXPECT_NEAR(slope[0], 2.0, 0.02);

	// l2-norm.
	EXPECT_NEAR(slope[1], 2.8, 0.03);

}


TEST(SecondOrderDerivative, EXPc4b4) {

	std::vector<double> slope = test_fd_approximation(0, 2, "d2dx2::c4b4", true, false);

	// Max norm.
	EXPECT_NEAR(slope[0], 2.0, 0.02);

	// l2-norm.
	EXPECT_NEAR(slope[1], 2.8, 0.03);

}


TEST(SecondOrderDerivative, COSc4b4) {

	std::vector<double> slope = test_fd_approximation(1, 2, "d2dx2::c4b4", true, false);

	// Max norm.
	EXPECT_NEAR(slope[0], 2.0, 0.02);

	// l2-norm.
	EXPECT_NEAR(slope[1], 2.8, 0.03);

}


TEST(SecondOrderDerivative, SUMc4b4) {

	std::vector<double> slope = test_fd_approximation(2, 2, "d2dx2::c4b4", true, false);

	// Max norm.
	EXPECT_NEAR(slope[0], 2.0, 0.02);

	// l2-norm.
	EXPECT_NEAR(slope[1], 2.8, 0.03);

}
