#include <cmath>


double y_func(
	const double kappa,
	const double sigma,
	const double time) {

	return sigma * sigma * (1 - std::exp(-2 * kappa * time)) / (2 * kappa);

}


double g_func(
	const double kappa,
	const double t_initial,
	const double t_final) {

	return (1 - std::exp(-kappa * (t_final - t_initial))) / kappa;

}
