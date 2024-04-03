#include <cmath>

#include "VasicekUtility.h"


double vasicek::b_func(
	const double time_1,
	const double time_2,
	const double kappa) {

	return (1 - std::exp(-kappa * (time_2 - time_1))) / kappa;

}


double vasicek::dbdt_func(
	const double time_1,
	const double time_2,
	const double kappa) {

	return -std::exp(-kappa * (time_2 - time_1));

}


double vasicek::a_func(
	const double time_1,
	const double time_2,
	const double kappa,
	const double theta,
	const double sigma) {

	const double sigma_sq = sigma * sigma;

	const double kappa_sq = kappa * kappa;

	const double b = vasicek::b_func(time_1, time_2, kappa);

	return (theta - sigma_sq / (2 * kappa_sq)) * (b - (time_2 - time_1)) 
		- sigma_sq * b * b / (4 * kappa);

}


double vasicek::dadt_func(
	const double time_1,
	const double time_2,
	const double kappa,
	const double theta,
	const double sigma) {


	const double sigma_sq = sigma * sigma;

	const double kappa_sq = kappa * kappa;

	const double b = vasicek::b_func(time_1, time_2, kappa);

	const double dbdt = vasicek::dbdt_func(time_1, time_2, kappa);

	return (theta - sigma_sq / (2 * kappa_sq)) * (dbdt + 1)
		- 2 * sigma_sq * b * dbdt / (4 * kappa);

}
