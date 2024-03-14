#pragma once

#include <complex>
#include <string>


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


std::complex<double> r_pm(
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
	const double lambda,
	const double theta,
	const double v_0,
	const double eta,
	const double rho,
	const double tau);
