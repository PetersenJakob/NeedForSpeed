#include "pch.h"


TEST(TestCaseName, TestName) {
  EXPECT_EQ(1, 1);
  EXPECT_TRUE(true);
}


TEST(FiniteDifference1stOrder, C2B1) {

	const int order = 21;

	const double x_max = 0.2;
	const double x_min = -0.2;

	const double dx = (x_max - x_min) / (order - 1);

	//
	std::vector<double> grid(order);

	for (int n = 0; n != order; ++n)
		grid.at(n) = dx * n + x_min;


//	TriDiagonal d1_c2b1_exp = d1dx1::c2b1(order, dx);


	std::vector<double> function(order, 0.0);
	std::vector<double> antideriv(order, 0.0);
	std::vector<double> antiantideriv(order, 0.0);
	std::vector<double> deriv(order, 0.0);
	std::vector<double> derivderiv(order, 0.0);



	for (int i = 0; i != order; ++i) {
		function[i] = cos(grid[i]);
		antideriv[i] = sin(grid[i]);
		antiantideriv[i] = -function[i];
		deriv[i] = -sin(grid[i]);
		derivderiv[i] = -function[i];
	}


//	std::vector<double> grid_exp = grid_equidistant(-0.2, 0.2, order);
//	std::vector<double> grid_cos = grid_equidistant(-M_PI / 2.0, M_PI / 2.0, order);

//	const double dx_exp = grid_exp[1] - grid_exp[0];
//	const double dx_cos = grid_cos[1] - grid_cos[0];

	EXPECT_EQ(1, 1);
	EXPECT_TRUE(true);
}
