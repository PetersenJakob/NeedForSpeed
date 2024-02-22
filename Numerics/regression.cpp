#include <vector>

#include "regression.h"


// Simple linear regression.
std::vector<double> regression::slr(
	const std::vector<double>& x, 
	const std::vector<double>& y) {

	double x_mean = 0.0;
	double y_mean = 0.0;

	double xx = 0.0;
	double xy = 0.0;

	for (int i = 0; i != x.size(); ++i) {

		x_mean += x[i];
		y_mean += y[i];

	}

	x_mean /= x.size();
	y_mean /= y.size();

	for (int i = 0; i != x.size(); ++i) {

		xx += pow(x[i] - x_mean, 2);
		xy += (x[i] - x_mean) * (y[i] - y_mean);

	}

	double beta = xy / xx;
	double alpha = y_mean - beta * x_mean;

	std::vector<double> result{ beta, alpha };

	return result;

}
