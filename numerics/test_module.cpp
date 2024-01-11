#include <cmath>
#include <iomanip>
#include <iostream>

#include "band_diagonal_matrix.h"
#include "finite_difference.h"
#include "grid_generator.h"


int main() {

	const int order = 21;

	TriDiagonal tri(order);
//	print_matrix(tri);

	PentaDiagonal penta(order);
//	print_matrix(penta);
	
	const double dx = 0.1;
	TriDiagonal d1dx1_ = d1dx1::c2(order, dx);
//	print_matrix(d1dx1_);

	d1dx1_.adjust_boundary();
	print_matrix(d1dx1_);

	TriDiagonal d2dx2_ = d2dx2::c2(order, dx);
//	print_matrix(d2dx2_);

	d2dx2_.adjust_boundary();
	print_matrix(d2dx2_);


	std::vector<double> grid = grid_equidistant(-1.0, 1.0, order);
	std::cout << std::scientific << std::setprecision(5);

	std::vector<double> func;
	for (int i = 0; i != order; ++i)
		func.push_back(exp(2 * grid[i]));

	std::vector<double> deriv1 = d1dx1_.mat_vec_product(func);

	for (int i = 0; i != order; ++i) {

		std::cout 
			<< std::setw(3) << i 
			<< std::setw(14) << grid[i] 
			<< std::setw(14) << func[i]
			<< std::setw(14) << 2 * exp(2 * grid[i])
			<< std::setw(14) << deriv1[i]
			<< std::endl;

	}

	return 0;
}