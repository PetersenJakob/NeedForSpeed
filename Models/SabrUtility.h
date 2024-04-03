#pragma once

#include <cmath>


namespace sabr {

	namespace implied_vol {

		double forward_mid(
			const double spot_forward,
			const double strike);

		double eta_func(
			const double spot_forward,
			const double spot_vol,
			const double strike,
			const double alpha,
			const double beta);

		double gamma_1(
			const double beta,
			const double f_mid);

		double gamma_2(
			const double beta,
			const double f_mid);

		double d_func(
			const double eta,
			const double rho);

		double black_scholes(
			const double spot_forward,
			const double spot_vol,
			const double alpha,
			const double beta,
			const double rho,
			const double strike,
			const double tau);

		double bachelier(
			const double spot_forward,
			const double spot_vol,
			const double alpha,
			const double beta,
			const double rho,
			const double strike,
			const double tau);

	}

	namespace call {

		double payoff(
			const double spot_forward,
			const double strike);

		double price(
			const double spot_forward,
			const double spot_vol,
			const double alpha,
			const double beta,
			const double rho,
			const double rate,
			const double strike,
			const double tau);

	}

}
