#include "pch.h"


TEST(CallOption, Test1) {

	double spot_price = 100.0;

	const double rate = 0.0;

	const double spot_variance = 0.01;

	const double lambda = 2.0;

	const double theta = 0.01;

	const double eta = 0.1;

	const double rho = 0.5;

	double strike = 100.0;

	const double tau = 0.5;

	double call_heston = 0.0;
	double call_bs = 0.0;

	double implied_vol = 0.0;

	std::cout << std::scientific << std::setprecision(6);

	for (int i = 0; i != 31; ++i) {

//		spot_price = 70.0 + 2 * i;

		strike = 70.0 + 2 * i;

		call_heston = heston::call(spot_price, spot_variance, rate, lambda, theta, eta, rho, strike, tau);

		call_bs = bs::call::price(spot_price, rate, std::sqrt(spot_variance), strike, tau);

//		implied_vol = bs::call::implied_vol(call_heston, spot_price, rate, strike, tau);

		std::cout
			<< std::setw(3) << i
			<< std::setw(16) << spot_price
			<< std::setw(16) << strike
			<< std::setw(16) << call_heston
			<< std::setw(16) << call_bs
			<< std::setw(16) << call_heston - call_bs
			<< std::setw(16) << implied_vol << std::endl;

	}

	std::cout << std::endl;

}
