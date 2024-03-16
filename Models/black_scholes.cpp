#include <cmath>

#include "black_scholes.h"
#include "distributions.h"


namespace bs {

	// European call option price.
	double call(
		const double spot_price,
		const double rate,
		const double sigma,
		const double strike,
		const double tau) {

		const double d_p = d_plus(spot_price, rate, sigma, strike, tau);

		const double d_m = d_minus(spot_price, rate, sigma, strike, tau);

		return normal::cdf(d_p) * spot_price 
			- normal::cdf(d_m) * strike * std::exp(-rate * tau);

	}

	// European put option price.
	double put(
		const double spot_price,
		const double rate,
		const double sigma,
		const double strike,
		const double tau) {

		const double call_price = call(spot_price, rate, sigma, strike, tau);
	
		// Put-call parity.
		return call_price - spot_price + strike * std::exp(-rate * tau);
	
	}

	namespace implied_vol {

		double call(
			const double spot_price,
			const double rate,
			const double sigma,
			const double strike,
			const double tau) {

			return 0.0;

		}

	}

}


double d_plus(
	const double spot_price,
	const double rate,
	const double sigma,
	const double strike,
	const double tau) {

	const double tmp = (rate + sigma * sigma / 2.0) * tau;

	return (std::log(spot_price / strike) + tmp) / (sigma * std::sqrt(tau));

}


double d_minus(
	const double spot_price,
	const double rate,
	const double sigma,
	const double strike,
	const double tau) {

	return d_plus(spot_price, rate, sigma, strike, tau) - sigma * std::sqrt(tau);

}
