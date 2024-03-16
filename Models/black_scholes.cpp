#include <cmath>

#include "black_scholes.h"
#include "distributions.h"


double bs::d_plus(
	const double spot_price,
	const double rate,
	const double sigma,
	const double strike,
	const double tau) {

	const double tmp = (rate + sigma * sigma / 2.0) * tau;

	return (std::log(spot_price / strike) + tmp) / (sigma * std::sqrt(tau));

}


double bs::d_minus(
	const double spot_price,
	const double rate,
	const double sigma,
	const double strike,
	const double tau) {

	return d_plus(spot_price, rate, sigma, strike, tau) - sigma * std::sqrt(tau);

}


// European call option price.
double bs::call::price(
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


// European call option delta.
double bs::call::delta(
	const double spot_price,
	const double rate,
	const double sigma,
	const double strike,
	const double tau) {

	const double d_p = d_plus(spot_price, rate, sigma, strike, tau);

	return normal::cdf(d_p);

}


// European call option gamma.
double bs::call::gamma(
	const double spot_price,
	const double rate,
	const double sigma,
	const double strike,
	const double tau) {

	const double d_p = d_plus(spot_price, rate, sigma, strike, tau);

	return normal::pdf(d_p) / (spot_price * sigma * std::sqrt(tau));

}


// European call option vega.
double bs::call::vega(
	const double spot_price,
	const double rate,
	const double sigma,
	const double strike,
	const double tau) {

	const double d_p = d_plus(spot_price, rate, sigma, strike, tau);

	return normal::pdf(d_p) * spot_price * std::sqrt(tau);

}


// European call option theta.
double bs::call::theta(
	const double spot_price,
	const double rate,
	const double sigma,
	const double strike,
	const double tau) {

	const double d_p = d_plus(spot_price, rate, sigma, strike, tau);

	const double d_m = d_minus(spot_price, rate, sigma, strike, tau);

	return -normal::pdf(d_p) * spot_price * sigma / (2.0 * std::sqrt(tau)) 
		- normal::cdf(d_m) * rate * strike * std::exp(-rate * tau);

}


// European call option rho.
double bs::call::rho(
	const double spot_price,
	const double rate,
	const double sigma,
	const double strike,
	const double tau) {

	const double d_m = d_minus(spot_price, rate, sigma, strike, tau);

	return normal::cdf(d_m) * strike * tau * std::exp(-rate * tau);

}


// Implied volatility corresponding to price of European call option.
double bs::call::implied_vol(
	const double option_price,
	const double spot_price,
	const double rate,
	const double strike,
	const double tau) {

	double sigma_1 = 0.0;
	double sigma_2 = 0.5;

	// Newton-Raphson root-search method.
	while (abs(sigma_2 - sigma_1) > 1.0e-5) {

		sigma_1 = sigma_2;

		double price = bs::call::price(spot_price, rate, sigma_1, strike, tau);
		double delta = bs::call::delta(spot_price, rate, sigma_1, strike, tau);

		sigma_2 = sigma_1 - (price - option_price) / delta;

	}

	return sigma_2;

}


// European put option price.
double bs::put::price(
	const double spot_price,
	const double rate,
	const double sigma,
	const double strike,
	const double tau) {

	const double call_price = 
		bs::call::price(spot_price, rate, sigma, strike, tau);

	// Put-call parity.
	return call_price - spot_price + strike * std::exp(-rate * tau);

}


// European put option delta.
double bs::put::delta(
	const double spot_price,
	const double rate,
	const double sigma,
	const double strike,
	const double tau) {

	return bs::call::delta(spot_price, rate, sigma, strike, tau) - 1.0;

}


// European put option gamma.
double bs::put::gamma(
	const double spot_price,
	const double rate,
	const double sigma,
	const double strike,
	const double tau) {

	return bs::call::gamma(spot_price, rate, sigma, strike, tau);

}


// European put option vega.
double bs::put::vega(
	const double spot_price,
	const double rate,
	const double sigma,
	const double strike,
	const double tau) {

	return bs::call::vega(spot_price, rate, sigma, strike, tau);

}


// European put option theta.
double bs::put::theta(
	const double spot_price,
	const double rate,
	const double sigma,
	const double strike,
	const double tau) {

	const double d_p = d_plus(spot_price, rate, sigma, strike, tau);

	const double d_m = d_minus(spot_price, rate, sigma, strike, tau);

	return -normal::pdf(d_p) * spot_price * sigma / (2.0 * std::sqrt(tau))
		+ normal::cdf(-d_m) * rate * strike * std::exp(-rate * tau);

}


// European put option rho.
double bs::put::rho(
	const double spot_price,
	const double rate,
	const double sigma,
	const double strike,
	const double tau) {

	const double d_m = d_minus(spot_price, rate, sigma, strike, tau);

	return -normal::cdf(-d_m) * strike * tau * std::exp(-rate * tau);

}


// Implied volatility corresponding to price of European put option.
double bs::put::implied_vol(
	const double option_price,
	const double spot_price,
	const double rate,
	const double strike,
	const double tau) {

	double sigma_1 = 0.0;
	double sigma_2 = 0.5;

	// Newton-Raphson root-search method.
	while (abs(sigma_2 - sigma_1) > 1.0e-5) {

		sigma_1 = sigma_2;

		double price = bs::put::price(spot_price, rate, sigma_1, strike, tau);
		double delta = bs::put::delta(spot_price, rate, sigma_1, strike, tau);

		sigma_2 = sigma_1 - (price - option_price) / delta;

	}

	return sigma_2;

}
