#pragma once


namespace normal {

	// Probability density function. Standard normal distribution.
	double pdf(const double x);

	// Probability density function. Normal distribution.
	double pdf(
		const double x,
		const double mu,
		const double sigma);

	// Cumulative distribution function. Standard normal distribution.
	double cdf(const double x);

	// Cumulative distribution function. Normal distribution.
	double cdf(
		const double x,
		const double mu,
		const double sigma);

}
