#pragma once

#include <complex>
#include <string>


namespace heston {

	double call(
		const double price,
		const double variance,
		const double rate,
		const double lambda,
		const double theta,
		const double eta,
		const double rho,
		const double strike,
		const double tau);

	double put(
		const double price,
		const double variance,
		const double rate,
		const double lambda,
		const double theta,
		const double eta,
		const double rho,
		const double strike,
		const double tau);

}


std::complex<double> alpha(
	const double j,
	const double k);


std::complex<double> beta(
	const double j,
	const double k,
	const double lambda,
	const double eta,
	const double rho);


double gamma(const double eta);


std::complex<double> discriminant(
	const double j,
	const double k,
	const double lambda,
	const double eta,
	const double rho);


std::complex<double> r_func(
	const std::string sign,
	const double j,
	const double k,
	const double lambda,
	const double eta,
	const double rho);


std::complex<double> g_func(
	const double j,
	const double k,
	const double lambda,
	const double eta,
	const double rho);


std::complex<double> d_func(
	const double j,
	const double k,
	const double lambda,
	const double eta,
	const double rho,
	const double tau);


std::complex<double> c_func(
	const double j,
	const double k,
	const double lambda,
	const double eta,
	const double rho,
	const double tau);


double probability(
	const double j,
	const double x,
	const double variance,
	const double lambda,
	const double theta,
	const double eta,
	const double rho,
	const double tau);
