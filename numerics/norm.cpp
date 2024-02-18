#include <algorithm>
#include <numeric>

#include "norm.h"


// Average function values...
std::vector<double> func_average(const std::vector<double>& func) {

	std::vector<double> result(func.size() - 1, 0.0);

	for (int i = 0; i != func.size() - 1; ++i) {
		result[i] = (func[i + 1] + func[i]) / 2.0;
	}

	return result;

}


// Average function values...
std::vector<std::vector<double>> 
func_average(const std::vector<std::vector<double>>& func) {

	std::vector<double> inner(func[0].size() - 1, 0.0);
	std::vector<std::vector<double>> result(func.size() - 1, inner);

	for (int i = 0; i != func.size() - 1; ++i) {
		for (int j = 0; j != func[0].size() - 1; ++j) {
			result[i][j] = 
				(func[i][j] + func[i + 1][j] 
					+ func[i][j + 1] + func[i + 1][j + 1]) / 4.0;
		}
	}

	return result;

}


// Infinity norm.
double norm::vector::infinity(const std::vector<double>& vec) {

	double norm = 0.0;

	for (int i = 0; i != vec.size(); ++i) {
		if (norm < abs(vec[i])) {
			norm = abs(vec[i]);
		}
	}

	return norm;

}


// l1 vector norm.
double norm::vector::l1(const std::vector<double>& vec) {

	double norm = 0.0;

	for (int i = 0; i != vec.size(); ++i) {
		norm += abs(vec[i]);
	}

	return norm;

}


// l2 vector norm.
double norm::vector::l2(const std::vector<double>& vec) {

	double norm = 0.0;

	for (int i = 0; i != vec.size(); ++i) {
		norm += vec[i] * vec[i];
	}

	return sqrt(norm);

}


// Infinity norm (1-dimension).
double norm::function::infinity(const std::vector<double>& func) {

	return norm::vector::infinity(func);

}


// Infinity norm (2-dimension).
double norm::function::infinity(const std::vector<std::vector<double>>& func) {

	double norm = 0.0;
	double norm_row = 0.0;

	for (int i = 0; i != func.size(); ++i) {
		norm_row = norm::function::infinity(func[i]);
		if (norm < norm_row) {
			norm = norm_row;
		}
	}

	return norm;

}


// L1 function norm (1-dimension).
double norm::function::l1(
	const double dx, 
	const std::vector<double>& func) {

	double norm = 0.0;

	std::vector<double> average = func_average(func);

	// Trapzoidal integration.
	for (int i = 0; i != func.size() - 1; ++i) {
		norm += dx * abs(average[i]);
	}

	return norm;

}


// L1 function norm (1-dimension).
double norm::function::l1(
	const std::vector<double>& grid, 
	const std::vector<double>& func) {

	double norm = 0.0;
	double dx = 0.0;

	std::vector<double> average = func_average(func);

	// Trapzoidal integration.
	for (int i = 0; i != func.size() - 1; ++i) {
		dx = grid[i + 1] - grid[i];
		norm += dx * abs(average[i]);
	}

	return norm;

}


// L1 function norm (2-dimension).
double norm::function::l1(
	const double dx,
	const double dy,
	const std::vector<std::vector<double>>& func) {

	double norm = 0.0;

	std::vector<std::vector<double>> average = func_average(func);

	// Trapzoidal integration.
	for (int i = 0; i != func.size() - 1; ++i) {
		for (int j = 0; j != func[0].size() - 1; ++j) {
			norm += dx * dy * abs(average[i][j]);
		}
	}

	return norm;

}


// L1 function norm (2-dimension).
double norm::function::l1(
	const std::vector<double>& grid_x,
	const std::vector<double>& grid_y,
	const std::vector<std::vector<double>>& func) {

	double norm = 0.0;
	double dx = 0.0;
	double dy = 0.0;

	std::vector<std::vector<double>> average = func_average(func);

	// Trapzoidal integration.
	for (int i = 0; i != func.size() - 1; ++i) {
		for (int j = 0; j != func[0].size() - 1; ++j) {
			dx = grid_x[i + 1] - grid_x[i];
			dy = grid_y[j + 1] - grid_y[j];
			norm += dx * dy * abs(average[i][j]);
		}
	}

	return norm;

}


// L2 function norm (1-dimension).
double norm::function::l2(
	const double dx, 
	const std::vector<double>& func) {

	double norm = 0.0;

	std::vector<double> average = func_average(func);

	// Trapzoidal integration.
	for (int i = 0; i != func.size() - 1; ++i) {
		norm += dx * pow(average[i], 2);
	}

	return sqrt(norm);

}


// L2 function norm (1-dimension).
double norm::function::l2(
	const std::vector<double>& grid, 
	const std::vector<double>& func) {

	double norm = 0.0;

	std::vector<double> average = func_average(func);

	// Trapzoidal integration.
	for (int i = 0; i != func.size() - 1; ++i) {
		norm += (grid[i + 1] - grid[i]) * pow(average[i], 2);
	}

	return sqrt(norm);

}


// L2 function norm (2-dimension).
double norm::function::l2(
	const double dx,
	const double dy,
	const std::vector<std::vector<double>>& func) {

	double norm = 0.0;

	std::vector<std::vector<double>> average = func_average(func);

	// Trapzoidal integration.
	for (int i = 0; i != func.size() - 1; ++i) {
		for (int j = 0; j != func[0].size() - 1; ++j) {
			norm += dx * dy * average[i][j] * average[i][j];
		}
	}

	return sqrt(norm);

}


// L2 function norm (2-dimension).
double norm::function::l2(
	const std::vector<double>& grid_x,
	const std::vector<double>& grid_y,
	const std::vector<std::vector<double>>& func) {

	double norm = 0.0;
	double dx = 0.0;
	double dy = 0.0;

	std::vector<std::vector<double>> average = func_average(func);

	// Trapzoidal integration.
	for (int i = 0; i != func.size() - 1; ++i) {
		for (int j = 0; j != func[0].size() - 1; ++j) {
			dx = grid_x[i + 1] - grid_x[i];
			dy = grid_y[j + 1] - grid_y[j];
			norm += dx * dy * average[i][j] * average[i][j];
		}
	}

	return sqrt(norm);

}


// Element-wise subtraction of vectors.
std::vector<double> norm::vector_diff(
	const std::vector<double>& vec1,
	const std::vector<double>& vec2) {

	std::vector<double> diff(vec1.size(), 0.0);

	for (int i = 0; i != vec1.size(); ++i) {
		diff[i] = vec1[i] - vec2[i];
	}

	return diff;

}


// Element-wise subtraction of matrices.
std::vector<std::vector<double>> norm::matrix_diff(
	const std::vector<std::vector<double>>& mat1,
	const std::vector<std::vector<double>>& mat2) {

	std::vector<double> inner(mat1[0].size(), 0.0);
	std::vector<std::vector<double>> diff(mat1.size(), inner);

	for (int i = 0; i != mat1.size(); ++i) {
		for (int j = 0; j != mat1[0].size(); ++j) {
			diff[i][j] = mat1[i][j] - mat2[i][j];
		}
	}

	return diff;

}
