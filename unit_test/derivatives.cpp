#include "pch.h"


std::vector<std::vector<double>> test_fd_approximation(
	const int function_type,
	const int derivative_order,
	const std::string fd_deriv_type,
	bool show_output_basic = false,
	bool show_output_all = false) {
	
	std::vector<double> dx_vec;
	std::vector<double> max_norm;
	std::vector<double> l2_norm;

	for (int i = 0; i != 15; ++i) {

		// Number of grid points.
		const int n_points = 21 + 10 * i;

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

	if (show_output_basic) {

		std::cout << std::scientific << std::setprecision(5);

		std::cout << std::endl << "               dx      max-norm       l2-norm" << std::endl;
		for (int i = 0; i != max_norm.size(); ++i) {

			std::cout
				<< std::setw(3) << i
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

	return ratio;

}


TEST(FirstOrderDerivative, EXPc2b1) {

	std::vector<std::vector<double>> ratio = test_fd_approximation(0, 1, "d1dx1::c2b1", false, false);

	for (int i = 0; i != ratio[0].size(); ++i) {
		EXPECT_NEAR(ratio[0][i], 2.0, 0.02);
	}

	for (int i = 0; i != ratio[1].size(); ++i) {
		EXPECT_NEAR(ratio[1][i], 2.8, 0.03);
	}

}


TEST(FirstOrderDerivative, COSc2b1) {

	std::vector<std::vector<double>> ratio = test_fd_approximation(1, 1, "d1dx1::c2b1", false, false);

	for (int i = 0; i != ratio[0].size(); ++i) {
		EXPECT_NEAR(ratio[0][i], 2.0, 0.02);
	}

	for (int i = 0; i != ratio[1].size(); ++i) {
		EXPECT_NEAR(ratio[1][i], 2.8, 0.05);
	}

}


TEST(FirstOrderDerivative, SUMc2b1) {

	std::vector<std::vector<double>> ratio = test_fd_approximation(2, 1, "d1dx1::c2b1", false, false);

	for (int i = 0; i != ratio[0].size(); ++i) {
		EXPECT_NEAR(ratio[0][i], 2.0, 0.02);
	}

	for (int i = 0; i != ratio[1].size(); ++i) {
		EXPECT_NEAR(ratio[1][i], 2.8, 0.07);
	}

}


TEST(FirstOrderDerivative, EXPc2b2) {

	std::vector<std::vector<double>> ratio = test_fd_approximation(0, 1, "d1dx1::c2b2", false, false);

	for (int i = 0; i != ratio[0].size(); ++i) {
		EXPECT_NEAR(ratio[0][i], 4.0, 0.06);
	}

	for (int i = 0; i != ratio[1].size(); ++i) {
		EXPECT_NEAR(ratio[1][i], 4.1, 0.3);
	}

}


TEST(FirstOrderDerivative, COSc2b2) {

	std::vector<std::vector<double>> ratio = test_fd_approximation(1, 1, "d1dx1::c2b2", false, false);

	for (int i = 0; i != ratio[0].size(); ++i) {
		EXPECT_NEAR(ratio[0][i], 3.9, 0.07);
	}

	for (int i = 0; i != ratio[1].size(); ++i) {
		EXPECT_NEAR(ratio[1][i], 4.3, 0.3);
	}

}


TEST(FirstOrderDerivative, SUMc2b2) {

	std::vector<std::vector<double>> ratio = test_fd_approximation(2, 1, "d1dx1::c2b2", false, false);

	for (int i = 0; i != ratio[0].size(); ++i) {
		EXPECT_NEAR(ratio[0][i], 3.9, 0.08);
	}

	for (int i = 0; i != ratio[1].size(); ++i) {
		EXPECT_NEAR(ratio[1][i], 4.3, 0.2);
	}

}


TEST(FirstOrderDerivative, EXPc4b4) {

	std::vector<std::vector<double>> ratio = test_fd_approximation(0, 1, "d1dx1::c4b4", false, false);

	for (int i = 0; i != ratio[0].size(); ++i) {
		EXPECT_NEAR(ratio[0][i], 16.0, 0.6);
	}

	for (int i = 0; i != ratio[1].size(); ++i) {
		EXPECT_NEAR(ratio[1][i], 20.0, 1.0);
	}

}


TEST(FirstOrderDerivative, COSc4b4) {

	std::vector<std::vector<double>> ratio = test_fd_approximation(1, 1, "d1dx1::c4b4", false, false);

	for (int i = 0; i != ratio[0].size(); ++i) {
		EXPECT_NEAR(ratio[0][i], 15.0, 0.7);
	}

	for (int i = 0; i != ratio[1].size(); ++i) {
		EXPECT_NEAR(ratio[1][i], 20.0, 0.6);
	}

}


TEST(FirstOrderDerivative, SUMc4b4) {

	std::vector<std::vector<double>> ratio = test_fd_approximation(2, 1, "d1dx1::c4b4", false, false);

	for (int i = 0; i != ratio[0].size(); ++i) {
		EXPECT_NEAR(ratio[0][i], 15.0, 0.8);
	}

	for (int i = 0; i != ratio[1].size(); ++i) {
		EXPECT_NEAR(ratio[1][i], 20.0, 0.6);
	}

}


TEST(SecondOrderDerivative, EXPc2b1) {

	std::vector<std::vector<double>> ratio = test_fd_approximation(0, 2, "d2dx2::c2b1", false, false);

	for (int i = 0; i != ratio[0].size(); ++i) {
		EXPECT_NEAR(ratio[0][i], 2.0, 0.03);
	}

	for (int i = 0; i != ratio[1].size(); ++i) {
		EXPECT_NEAR(ratio[1][i], 2.8, 0.03);
	}

}


TEST(SecondOrderDerivative, COSc2b1) {

	std::vector<std::vector<double>> ratio = test_fd_approximation(1, 2, "d2dx2::c2b1", false, false);

	for (int i = 0; i != ratio[0].size(); ++i) {
		EXPECT_NEAR(ratio[0][i], 2.0, 0.06);
	}

	for (int i = 0; i != ratio[1].size(); ++i) {
		EXPECT_NEAR(ratio[1][i], 2.8, 0.05);
	}

}


TEST(SecondOrderDerivative, SUMc2b1) {

	std::vector<std::vector<double>> ratio = test_fd_approximation(2, 2, "d2dx2::c2b1", false, false);

	for (int i = 0; i != ratio[0].size(); ++i) {
		EXPECT_NEAR(ratio[0][i], 2.0, 0.05);
	}

	for (int i = 0; i != ratio[1].size(); ++i) {
		EXPECT_NEAR(ratio[1][i], 2.8, 0.04);
	}

}


TEST(SecondOrderDerivative, EXPc2b2) {

	std::vector<std::vector<double>> ratio = test_fd_approximation(0, 2, "d2dx2::c2b2", false, false);

	for (int i = 0; i != ratio[0].size(); ++i) {
		EXPECT_NEAR(ratio[0][i], 4.0, 0.09);
	}

	for (int i = 0; i != ratio[1].size(); ++i) {
		EXPECT_NEAR(ratio[1][i], 5.3, 0.2);
	}

}


TEST(SecondOrderDerivative, COSc2b2) {

	std::vector<std::vector<double>> ratio = test_fd_approximation(1, 2, "d2dx2::c2b2", false, false);

	for (int i = 0; i != ratio[0].size(); ++i) {
		EXPECT_NEAR(ratio[0][i], 4.0, 0.09);
	}

	for (int i = 0; i != ratio[1].size(); ++i) {
		EXPECT_NEAR(ratio[1][i], 5.25, 0.3);
	}

}


TEST(SecondOrderDerivative, SUMc2b2) {

	std::vector<std::vector<double>> ratio = test_fd_approximation(2, 2, "d2dx2::c2b2", false, false);

	for (int i = 0; i != ratio[0].size(); ++i) {
		EXPECT_NEAR(ratio[0][i], 4.0, 0.05);
	}

	for (int i = 0; i != ratio[1].size(); ++i) {
		EXPECT_NEAR(ratio[1][i], 5.25, 0.3);
	}

}


TEST(SecondOrderDerivative, EXPc4b4) {

	std::vector<std::vector<double>> ratio = test_fd_approximation(0, 2, "d2dx2::c4b4", false, false);

	for (int i = 0; i != ratio[0].size(); ++i) {
		EXPECT_NEAR(ratio[0][i], 15.0, 2.0);
	}

	for (int i = 0; i != ratio[1].size(); ++i) {
		EXPECT_NEAR(ratio[1][i], 20.0, 2.5);
	}

}


TEST(SecondOrderDerivative, COSc4b4) {

	std::vector<std::vector<double>> ratio = test_fd_approximation(1, 2, "d2dx2::c4b4", false, false);

	for (int i = 0; i != ratio[0].size(); ++i) {
		EXPECT_NEAR(ratio[0][i], 16.0, 0.9);
	}

	for (int i = 0; i != ratio[1].size(); ++i) {
		EXPECT_NEAR(ratio[1][i], 23.0, 0.7);
	}

}


TEST(SecondOrderDerivative, SUMc4b4) {

	std::vector<std::vector<double>> ratio = test_fd_approximation(2, 2, "d2dx2::c4b4", false, false);

	for (int i = 0; i != ratio[0].size(); ++i) {
		EXPECT_NEAR(ratio[0][i], 16.0, 0.9);
	}

	for (int i = 0; i != ratio[1].size(); ++i) {
		EXPECT_NEAR(ratio[1][i], 23.0, 0.8);
	}

}
