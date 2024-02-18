#pragma once

#include <vector>


namespace norm {

	namespace vector {

		// Infinity norm.
		double infinity(std::vector<double> vec);

		// l1 vector norm.
		double l1(std::vector<double> vec);

		// l2 vector norm.
		double l2(std::vector<double> vec);

	}

	namespace function {

		// Infinity norm (1-dimension).
		double infinity(std::vector<double> func);

		// Infinity norm (2-dimension).
		double infinity(std::vector<std::vector<double>> func);

		// L1 function norm (1-dimension).
		double l1(
			const double dx, 
			std::vector<double> func);

		// L1 function norm (1-dimension).
		double l1(
			const std::vector<double>& grid, 
			std::vector<double> func);

		// L1 function norm (2-dimension).
		double l1(
			const double dx,
			const double dy,
			std::vector<std::vector<double>> func);

		// L1 function norm (2-dimension).
		double l1(
			const std::vector<double>& grid_x,
			const std::vector<double>& grid_y,
			std::vector<std::vector<double>> func);

		// L2 function norm (1-dimension).
		double l2(
			const double dx, 
			std::vector<double> func);

		// L2 function norm (1-dimension).
		double l2(
			const std::vector<double>& grid,
			std::vector<double> func);

		// L2 function norm (2-dimension).
		double l2(
			const double dx,
			const double dy,
			std::vector<std::vector<double>> func);

		// L2 function norm (2-dimension).
		double l2(
			const std::vector<double>& grid_x,
			const std::vector<double>& grid_y,
			std::vector<std::vector<double>> func);

	}

	// Element-wise subtraction of vectors.
	std::vector<double> vector_diff(
		const std::vector<double>& vec1,
		const std::vector<double>& vec2);

	// Element-wise subtraction of matrices.
	std::vector<std::vector<double>> matrix_diff(
		const std::vector<std::vector<double>>& mat1,
		const std::vector<std::vector<double>>& mat2);

}
