#pragma once


namespace bs {

	// European call option price.
	double call(
		const double spot_price,
		const double rate,
		const double sigma,
		const double strike,
		const double tau);

	// European put option price.
	double put(
		const double spot_price,
		const double rate,
		const double sigma,
		const double strike,
		const double tau);

	namespace implied_vol {

		double call(
			const double spot_price,
			const double rate,
			const double sigma,
			const double strike,
			const double tau);

	}

}


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
