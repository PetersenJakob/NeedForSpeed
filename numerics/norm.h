#pragma once

#include <vector>


namespace norm {

	namespace vector {

		// Maximum norm.
		double max(std::vector<double> vec);

		// l1 vector norm.
		double l1(std::vector<double> vec);

		// l2 vector norm.
		double l2(std::vector<double> vec);

	}

	namespace function {

		// L1 function norm. 1-dimensional.
		double l1(
			const double dx, 
			std::vector<double> vec);

		// L1 function norm. 1-dimensional.
		double l1(
			const std::vector<double>& grid, 
			std::vector<double> vec);

		// L1 function norm. 2-dimensional.
		double l1(
			const double dx,
			const double dy,
			std::vector<std::vector<double>> vec);

		// L1 function norm. 2-dimensional.
		double l1(
			const std::vector<double>& grid_x,
			const std::vector<double>& grid_y,
			std::vector<std::vector<double>> vec);

		// L2 function norm. 1-dimensional.
		double l2(
			const double dx, 
			std::vector<double> vec);

		// L2 function norm. 1-dimensional.
		double l2(
			const std::vector<double>& grid,
			std::vector<double> vec);

		// L2 function norm. 2-dimensional.
		double l2(
			const double dx,
			const double dy,
			std::vector<std::vector<double>> vec);

		// L2 function norm. 2-dimensional.
		double l2(
			const std::vector<double>& grid_x,
			const std::vector<double>& grid_y,
			std::vector<std::vector<double>> vec);

	}

	// Element-wise subtraction of vectors.
	std::vector<double> vector_diff(
		const std::vector<double>& vec1,
		const std::vector<double>& vec2);

	// Element-wise subtraction of matrices.
	std::vector<std::vector<double>> vector_diff(
		const std::vector<std::vector<double>>& mat1,
		const std::vector<std::vector<double>>& mat2);

}
