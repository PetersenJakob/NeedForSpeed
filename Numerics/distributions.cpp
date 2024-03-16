#define _USE_MATH_DEFINES
#include <cmath>

#include "distributions.h"


namespace normal {

	// Probability density function. Standard normal distribution.
	double pdf(const double x) {

		return std::exp(-x * x / 2.0) / std::sqrt(2.0 * M_PI);

	}

	// Probability density function. Normal distribution.
	double pdf(
		const double x,
		const double mu,
		const double sigma) {

		const double x_tmp = (x - mu) / sigma;

		return std::exp(-x_tmp * x_tmp / 2.0) / (sigma * std::sqrt(2.0 * M_PI));

	}

	// Cumulative distribution function. Standard normal distribution.
	double cdf(const double x) {

		return std::erfc(-x / std::sqrt(2.0)) / 2.0;

	}

	// Cumulative distribution function. Normal distribution.
	double cdf(
		const double x,
		const double mu,
		const double sigma) {

		return std::erfc(-x / std::sqrt(2.0)) / 2.0;

	}

}
