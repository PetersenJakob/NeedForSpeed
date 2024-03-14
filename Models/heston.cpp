#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>

#include "heston.h"


// Heston model, Gatheral 2006
//	dS_t = rate * S_t * dt + \sqrt{v_t} * S_t * dW1_t
//	dv_t = - lambda * (v_t - theta) * dt + eta * \sqrt{v_t} * dW2_t
//  <dW1_t, dW2_t> = rho * dt


std::complex<double> alpha(
	const double j,
	const double k) {

	std::complex<double> i(0.0, 1.0);

	return - k * k / 2.0 - i * k / 2.0 + i * j * k;

}


std::complex<double> beta(
	const double j,
	const double k,
	const double lambda,
	const double eta,
	const double rho) {

	std::complex<double> i(0.0, 1.0);

	return lambda - rho * eta * j - i * rho * eta * k;

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


std::complex<double> r_pm(
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

	std::complex<double> r_minus = r_pm("minus", j, k, lambda, eta, rho);

	std::complex<double> r_plus = r_pm("plus", j, k, lambda, eta, rho);

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

	std::complex<double> r_minus = r_pm("minus", j, k, lambda, eta, rho);

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

	std::complex<double> r_minus = r_pm("minus", j, k, lambda, eta, rho);

	const double gamma_ = gamma(eta);

	std::complex<double> result = 1.0 - g * std::exp(-d * tau);

	result /= 1.0 - g;

	result = lambda * (r_minus * tau - log(result) / gamma);

	return result;

}


double probability(
	const double j,
	const double x,
	const double lambda,
	const double theta,
	const double v_0,
	const double eta,
	const double rho,
	const double tau) {

	std::complex<double> c(0.0, 0.0);
	std::complex<double> d(0.0, 0.0);
	
	const int n_steps = 100;

	const double k_max = 100.0;

	const double increment = k_max / (n_steps - 1);

	double integral = 0.0;

	for (int i = 0; i != n_steps; ++i) {

		double k = increment * (i + 0.5);

		c = c_func(j, k, lambda, eta, rho, tau);
		d = d_func(j, k, lambda, eta, rho, tau);

		integral += 
			std::real(std::exp(c * theta + d * v_0 + i * k * x) / (i * k));

	}

	return 0.5 + integral / M_PI;

}
