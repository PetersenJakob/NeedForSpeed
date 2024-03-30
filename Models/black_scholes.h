#pragma once

#include <functional>
#include <vector>




namespace bs {

	double d_plus(
		const double spot_price,
		const double rate,
		const double sigma,
		const double strike,
		const double tau);

	double d_minus(
		const double spot_price,
		const double rate,
		const double sigma,
		const double strike,
		const double tau);

	namespace pde {

		namespace generator {

			std::vector<std::vector<double>> prefactor(
				const double rate,
				const double sigma,
				const std::vector<double>& spatial_grid);

			template <class T>
			T derivative_full(
				const double rate,
				const double sigma,
				const std::vector<double>& spatial_grid,
				const std::vector<std::function<T(std::vector<double>)>>& deriv) {

				std::vector<std::vector<double>> prefactor_ 
					= prefactor(rate, sigma, spatial_grid);

				// Identity operator.
				T identity = deriv[0](spatial_grid).identity();
				T derivative = identity.pre_vector(prefactor_[0]);
				// First order derivative operator.
				derivative += deriv[0](spatial_grid).pre_vector(prefactor_[1]);
				// Second order derivative operator.
				derivative += deriv[1](spatial_grid).pre_vector(prefactor_[2]);

				return derivative;

			}

			template <class T>
			std::function<T(const std::vector<double>&)> derivative(
				const double rate,
				const double sigma,
				const std::vector<std::function<T(std::vector<double>)>>& deriv) {

				return [rate, sigma, deriv](
					const std::vector<double>& spatial_grid) {
						return derivative_full<T>(rate, sigma, spatial_grid, deriv);
					};

			}

		}

	}

	// European call option.
	namespace call {

		double payoff(
			const double spot_price,
			const double strike);

		double price(
			const double spot_price,
			const double rate,
			const double sigma,
			const double strike,
			const double tau);

		double delta(
			const double spot_price,
			const double rate,
			const double sigma,
			const double strike,
			const double tau);

		double gamma(
			const double spot_price,
			const double rate,
			const double sigma,
			const double strike,
			const double tau);

		double vega(
			const double spot_price,
			const double rate,
			const double sigma,
			const double strike,
			const double tau);

		double theta(
			const double spot_price,
			const double rate,
			const double sigma,
			const double strike,
			const double tau);

		double rho(
			const double spot_price,
			const double rate,
			const double sigma,
			const double strike,
			const double tau);

		double implied_vol(
			const double option_price,
			const double spot_price,
			const double rate,
			const double strike,
			const double tau);

		std::function<std::vector<double>
			(const double, const std::vector<std::vector<double>>&)>
			solution_func(
				const double rate,
				const double sigma,
				const double strike);

		std::vector<double> solution_full(
			const std::vector<double>& spatial_grid,
			const double rate,
			const double sigma,
			const double strike,
			const double tau);

	}

	// European put option.
	namespace put {

		double price(
			const double spot_price,
			const double rate,
			const double sigma,
			const double strike,
			const double tau);

		double delta(
			const double spot_price,
			const double rate,
			const double sigma,
			const double strike,
			const double tau);

		double gamma(
			const double spot_price,
			const double rate,
			const double sigma,
			const double strike,
			const double tau);

		double vega(
			const double spot_price,
			const double rate,
			const double sigma,
			const double strike,
			const double tau);

		double theta(
			const double spot_price,
			const double rate,
			const double sigma,
			const double strike,
			const double tau);

		double rho(
			const double spot_price,
			const double rate,
			const double sigma,
			const double strike,
			const double tau);

		double implied_vol(
			const double option_price,
			const double spot_price,
			const double rate,
			const double strike,
			const double tau);

	}

}
