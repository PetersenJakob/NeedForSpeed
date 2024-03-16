#include "pch.h"


TEST(CallOption, Test1) {

	const double spot_price = 100.0;

	const double rate = 0.02;

	const double spot_variance = 0.01;

	const double lambda = 2.0;

	const double theta = 0.01;

	const double eta = 0.1;

	const double rho = 0.5;

	double strike = 0.0;

	const double tau = 0.5;

	double call_price = 0.0;

	std::cout << std::scientific << std::setprecision(8);

	for (int i = 0; i != 31; ++i) {

		strike = 85.0 + i;

		call_price = heston::call(spot_price, spot_variance, rate, lambda, theta, eta, rho, strike, tau);

		std::cout
			<< std::setw(3) << i
			<< std::setw(20) << strike
			<< std::setw(20) << call_price << std::endl;

	}
	std::cout << std::endl;

}
