#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>

#include "heston.h"


// Heston model, see Gatheral (2006).
// The underlying price process is given by
//		dS_t = rate * S_t * dt + \sqrt{v_t} * S_t * dW1_t,
// where the variance process reads
//		dv_t = lambda * (theta - v_t) * dt + eta * \sqrt{v_t} * dW2_t
// and the Wiener processes are correlated
//		<dW1_t, dW2_t> = rho * dt.
//


namespace heston {

	double call(
		const double spot_price,
		const double spot_variance,
		const double rate,
		const double lambda,
		const double theta,
		const double eta,
		const double rho,
		const double strike,
		const double tau) {

		const double forward_price = spot_price * std::exp(rate * tau);

		const double x = std::log(forward_price / strike);

		const double prop_0 = 
			probability(0, x, spot_variance, lambda, theta, eta, rho, tau);

		const double prop_1 =
			probability(1, x, spot_variance, lambda, theta, eta, rho, tau);

		return std::exp(x) * prop_1 - prop_0;
	}

	double put(
		const double spot_price,
		const double spot_variance,
		const double rate,
		const double lambda,
		const double theta,
		const double eta,
		const double rho,
		const double strike,
		const double tau) {

		const double call_price = 
			call(spot_price, spot_variance, rate, lambda, theta, eta, rho, strike, tau);

		// Put-call parity.
		return call_price - spot_price + strike * std::exp(-rate * tau);

	}

}


std::complex<double> alpha(
	const double j,
	const double k) {

	std::complex<double> i_unit(0.0, 1.0);

	return - k * k / 2.0 - i_unit * k / 2.0 + i_unit * j * k;

}


std::complex<double> beta(
	const double j,
	const double k,
	const double lambda,
	const double eta,
	const double rho) {

	std::complex<double> i_unit(0.0, 1.0);

	return lambda - rho * eta * j - i_unit * rho * eta * k;

}


double gamma(const double eta) {

	return eta * eta / 2.0;

}


std::complex<double> discriminant(
	const double j,
	const double k,
	const double lambda,
	const double eta,
	const double rho) {

	const std::complex<double> alpha_ = alpha(j, k);

	const std::complex<double> beta_ = beta(j, k, lambda, eta, rho);

	const double gamma_ = gamma(eta);

	return std::sqrt(beta_ * beta_ - 4.0 * alpha_ * gamma_);

}


std::complex<double> r_func(
	const std::string sign,
	const double j,
	const double k,
	const double lambda,
	const double eta,
	const double rho) {

	const std::complex<double> beta_ = beta(j, k, lambda, eta, rho);

	const double gamma_ = gamma(eta);

	const std::complex<double> d = discriminant(j, k, lambda, eta, rho);

	if (sign == "plus") {
		return (beta_ + d) / (2.0 * gamma_);
	}
	else if (sign == "minus") {
		return (beta_ - d) / (2.0 * gamma_);
	}
	else {
		throw std::invalid_argument("sign unknown.");
	}

}


std::complex<double> g_func(
	const double j,
	const double k,
	const double lambda,
	const double eta,
	const double rho) {

	std::complex<double> r_minus = r_func("minus", j, k, lambda, eta, rho);

	std::complex<double> r_plus = r_func("plus", j, k, lambda, eta, rho);

	return r_minus / r_plus;

}


std::complex<double> d_func(
	const double j,
	const double k,
	const double lambda,
	const double eta,
	const double rho,
	const double tau) {

	const std::complex<double> d = discriminant(j, k, lambda, eta, rho);

	const std::complex<double> g = g_func(j, k, lambda, eta, rho);

	std::complex<double> r_minus = r_func("minus", j, k, lambda, eta, rho);

	return r_minus * (1.0 - std::exp(-d * tau)) / (1.0 - g * std::exp(-d * tau));

}


std::complex<double> c_func(
	const double j,
	const double k,
	const double lambda,
	const double eta,
	const double rho,
	const double tau) {

	const std::complex<double> d = discriminant(j, k, lambda, eta, rho);

	const std::complex<double> g = g_func(j, k, lambda, eta, rho);

	std::complex<double> r_minus = r_func("minus", j, k, lambda, eta, rho);

	const double gamma_ = gamma(eta);

	std::complex<double> result = 1.0 - g * std::exp(-d * tau);

	result /= 1.0 - g;

	return lambda * (r_minus * tau - std::log(result) / gamma_);

}


double probability(
	const double j,
	const double x,
	const double spot_variance,
	const double lambda,
	const double theta,
	const double eta,
	const double rho,
	const double tau) {

	// Number of integration steps.	
	const int n_steps = 101;

	// Integration range [0; k_max].
	const double k_max = 100.0;

	// Integration step size.
	const double step_size = k_max / (n_steps - 1);

	double integral = 0.0;
	std::complex<double> c(0.0, 0.0);
	std::complex<double> d(0.0, 0.0);

	std::complex<double> i_unit(0.0, 1.0);

	// Integral represented by simple Riemann sum.
	for (int i = 0; i != n_steps - 1; ++i) {

		double k = step_size * (i + 0.5);

		c = c_func(j, k, lambda, eta, rho, tau);
		d = d_func(j, k, lambda, eta, rho, tau);

		c = std::exp(c * theta + d * spot_variance + i_unit * k * x);
		c /= i_unit * k;

		integral += std::real(c) * step_size;

	}

	return 0.5 + integral / M_PI;

}
