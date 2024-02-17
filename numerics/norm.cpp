#include <algorithm>
#include <numeric>

#include "norm.h"


// Maximum norm.
double norm::vector::max(std::vector<double> vec) {

	// Absolue value of each element.
	std::transform(vec.begin(), vec.end(), vec.begin(), [](double x) { return abs(x); });

	return *std::max_element(vec.begin(), vec.end());

}


// l1 vector norm.
double norm::vector::l1(std::vector<double> vec) {

	// Absolute value of each element.
	std::transform(vec.begin(), vec.end(), vec.begin(), [](double x) { return abs(x); });

	return std::accumulate(vec.begin(), vec.end(), (double)0);

}


// l2 vector norm.
double norm::vector::l2(std::vector<double> vec) {

	// Square of each element.
	std::transform(vec.begin(), vec.end(), vec.begin(), [](double x) { return x * x; });

	return std::accumulate(vec.begin(), vec.end(), (double)0);

}


// L1 function norm. 1-dimensional.
double norm::function::l1(
	const double dx, 
	std::vector<double> vec) {

	double norm = 0.0;

	// Trapzoidal integration.
	for (int i = 0; i != vec.size() - 1; ++i) {
		norm += dx * abs((vec[i + 1] + vec[i]) / 2.0);
	}

	return norm;

}


// L1 function norm. 1-dimensional.
double norm::function::l1(
	const std::vector<double>& grid, 
	std::vector<double> vec) {

	double norm = 0.0;

	// Trapzoidal integration.
	for (int i = 0; i != vec.size() - 1; ++i) {
		norm += (grid[i + 1] - grid[i]) * abs((vec[i + 1] + vec[i]) / 2.0);
	}

	return norm;

}


// TODO: Should return the square root!
// L2 function norm. 1-dimensional.
double norm::function::l2(
	const double dx, 
	std::vector<double> vec) {

	double norm = 0.0;

	// Trapzoidal integration.
	for (int i = 0; i != vec.size() - 1; ++i) {
		norm += dx * pow((vec[i + 1] + vec[i]) / 2.0, 2);
	}

	return norm;

}


// TODO: Should return the square root!
// L2 function norm. 1-dimensional.
double norm::function::l2(
	const std::vector<double>& grid, 
	std::vector<double> vec) {

	double norm = 0.0;

	// Trapzoidal integration.
	for (int i = 0; i != vec.size() - 1; ++i) {
		norm += (grid[i + 1] - grid[i]) * pow((vec[i + 1] + vec[i]) / 2.0, 2);
	}

	return norm;

}


// Element-wise subtraction of vectors.
std::vector<double> norm::vector_diff(
	const std::vector<double>& vec1,
	const std::vector<double>& vec2) {

	std::vector<double> diff(vec1.size(), 0.0);

	std::transform(vec1.begin(), vec1.end(), vec2.begin(), diff.begin(), std::minus<double>());

	return diff;

}