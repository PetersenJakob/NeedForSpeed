#include "pch.h"


TEST(CallOption, Test1) {

	double spot_price = 100.0;

	const double rate = 0.0;

	const double spot_variance_1 = 0.01;
	const double spot_variance_2 = 0.011;
	const double spot_variance_3 = 0.01;

	const double lambda_1 = 2.0;
	const double lambda_2 = 4.0;
	const double lambda_3 = 2.0;

	const double theta = 0.01;

	const double eta_1 = 0.1;
	const double eta_2 = 0.08;
	const double eta_3 = 0.1;

	const double rho_1 =  0.0;
	const double rho_2 = -0.25;
	const double rho_3 =  0.0;

	double strike = 100.0;

	const double tau = 0.5;

	std::ofstream myfile("heston.csv");
	myfile << std::scientific << std::setprecision(6);

	for (int i = 0; i != 61; ++i) {

//		spot_price = 70.0 + 2 * i;

		strike = 70.0 + 1 * i;

		double call_heston_1 = heston::call(spot_price, spot_variance_1, rate, lambda_1, theta, eta_1, rho_1, strike, tau);
		double call_heston_2 = heston::call(spot_price, spot_variance_2, rate, lambda_2, theta, eta_2, rho_2, strike, tau);
		double call_heston_3 = heston::call(spot_price, spot_variance_3, rate, lambda_3, theta, eta_3, rho_3, strike, tau);

		double call_bs_1 = bs::call::price(spot_price, rate, std::sqrt(spot_variance_1), strike, tau);
		double call_bs_2 = bs::call::price(spot_price, rate, std::sqrt(spot_variance_2), strike, tau);
		double call_bs_3 = bs::call::price(spot_price, rate, std::sqrt(spot_variance_3), strike, tau);

		double implied_vol_1 = bs::call::implied_vol(call_heston_1, spot_price, rate, strike, tau);
		double implied_vol_2 = bs::call::implied_vol(call_heston_2, spot_price, rate, strike, tau);

		myfile
			<< std::setw(3) << i << ", "
			<< std::setw(16) << spot_price << ", "
			<< std::setw(16) << strike << ", "
			<< std::setw(16) << call_heston_1 - call_bs_1 << ", "
			<< std::setw(16) << call_heston_2 - call_bs_2 << ", "
			<< std::setw(16) << call_heston_3 - call_bs_3 << ", "
			<< std::setw(16) << implied_vol_1 << ", "
			<< std::setw(16) << implied_vol_2 << std::endl;

	}

	myfile.close();

}


TEST(CallOption, FD1) {

	const double s_min = 0.0;
	const double s_max = 200.0;

	const double v_min = 0.0;
	const double v_max = 1.0;

	const double rate = 0.03;

	const double strike = 100.0;

	const double tau = 1.0;

	const double lambda = 3.0;
	const double theta = 0.12;
	const double eta = 0.041;
	const double rho = 0.0;

	std::vector<double> dx_vec;
	std::vector<double> dy_vec;
	std::vector<double> max_norm;
	std::vector<double> l1_norm;
	std::vector<double> l2_norm;

	const int n_iterations = 1;

	for (int i = 0; i != n_iterations; ++i) {

		const int n_points_s = 51 + i * 5;
		const int n_points_v = 51 + i * 5;

		const std::vector<double> grid_s = grid::uniform(s_min, s_max, n_points_s);
		const std::vector<double> grid_v = grid::uniform(v_min, v_max, n_points_v);


		// S-dimension.
		TriDiagonal d1_s = d1dx1::uniform::c2b1(grid_s);
		d1_s.pre_vector(grid_s);
		d1_s *= rate;

		TriDiagonal d2_s = d2dx2::uniform::c2b0(grid_s);
		d2_s.pre_vector(grid_s);
		d2_s.pre_vector(grid_s);
		d2_s *= rate;



	}

}
