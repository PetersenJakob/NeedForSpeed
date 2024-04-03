#include <cmath>

#include "BlackScholesUtility.h"
#include "SabrUtility.h"


double sabr::implied_vol::forward_mid(
	const double spot_forward,
	const double strike) {

	return (spot_forward + strike) / 2;

}


double sabr::implied_vol::eta_func(
	const double spot_forward,
	const double spot_vol,
	const double strike,
	const double alpha,
	const double beta) {

	const double f_beta = std::pow(spot_forward, 1 - beta);
	const double k_beta = std::pow(strike, 1 - beta);

	return alpha * (f_beta - k_beta) / (spot_vol * (1 - beta));

}


double sabr::implied_vol::gamma_1(
	const double beta,
	const double f_mid) {

	return beta / f_mid;

}


double sabr::implied_vol::gamma_2(
	const double beta,
	const double f_mid) {

	return -beta * (1 - beta) / std::pow(f_mid, 2);

}


double sabr::implied_vol::d_func(
	const double eta,
	const double rho) {

	const double numerator = std::sqrt(1 - 2 * rho * eta + eta * eta);

	return std::log((numerator + eta - rho) / (1 - rho));

}


double sabr::implied_vol::black_scholes(
	const double spot_forward,
	const double spot_vol,
	const double alpha,
	const double beta,
	const double rho,
	const double strike,
	const double tau) {

	const double f_mid = sabr::implied_vol::forward_mid(spot_forward, strike);

	const double eta = 
		sabr::implied_vol::eta_func(spot_forward, spot_vol, strike, alpha, beta);

	const double g_1 = sabr::implied_vol::gamma_1(beta, f_mid);

	const double g_2 = sabr::implied_vol::gamma_2(beta, f_mid);

	const double d = sabr::implied_vol::d_func(eta, rho);

	double implied_vol = (2 * g_2 - g_1 * g_1 + 1 / (f_mid * f_mid)) / 24;
	implied_vol *= std::pow(spot_vol * std::pow(f_mid, beta) / alpha, 2);
	implied_vol += rho * g_1 * spot_vol * std::pow(f_mid, beta) / (4 * alpha);
	implied_vol += (2 - 3 * rho * rho) / 24;
	implied_vol *= tau * alpha * alpha;
	implied_vol += 1;
	implied_vol *= alpha * std::log(spot_forward / strike) / d;

	return implied_vol;

}


double sabr::implied_vol::bachelier(
	const double spot_forward,
	const double spot_vol,
	const double alpha,
	const double beta,
	const double rho,
	const double strike,
	const double tau) {

	const double f_mid = sabr::implied_vol::forward_mid(spot_forward, strike);

	const double eta = 
		sabr::implied_vol::eta_func(spot_forward, spot_vol, strike, alpha, beta);

	const double g_1 = sabr::implied_vol::gamma_1(beta, f_mid);

	const double g_2 = sabr::implied_vol::gamma_2(beta, f_mid);

	const double d = sabr::implied_vol::d_func(eta, rho);

	double implied_vol = (2 * g_2 - g_1 * g_1) / 24;
	implied_vol *= std::pow(spot_vol * std::pow(f_mid, beta) / alpha, 2);
	implied_vol += rho * g_1 * spot_vol * std::pow(f_mid, beta) / (4 * alpha);
	implied_vol += (2 - 3 * rho * rho) / 24;
	implied_vol *= tau * alpha * alpha;
	implied_vol += 1;
	implied_vol *= alpha * (spot_forward - strike) / d;

	return implied_vol;

}


double sabr::call::payoff(
	const double spot_forward,
	const double strike) {

	bs::call::payoff(spot_forward, strike);

}


double sabr::call::price(
	const double spot_forward,
	const double spot_vol,
	const double alpha,
	const double beta,
	const double rho,
	const double rate,
	const double strike,
	const double tau) {

	const double implied_vol = 
		sabr::implied_vol::black_scholes(spot_forward, spot_vol, alpha, beta, rho, strike, tau);

	return bs::call::price(spot_forward, rate, implied_vol, strike, tau);

}
