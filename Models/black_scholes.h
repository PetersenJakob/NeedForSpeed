#pragma once


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

	// European call option.
	namespace call {

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
