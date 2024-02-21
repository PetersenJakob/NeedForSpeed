#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "derivatives.h"

int main() {

	const double x_min = -1.0;
	const int size = 21;
	const double dx = 0.1;

	std::vector<element> deriv_operator = d1dx1::central_2o_3p(size, dx);

	std::cout << std::scientific << std::setprecision(6);

	std::vector<double> grid;
	grid.reserve(size);
	for (int n = size - 1; n >= 0; --n)
		grid.push_back(n * dx + x_min);

	std::vector<double> func;
	func.reserve(size);
	std::vector<double> deriv;
	deriv.reserve(size);
	for (int n = 0; n != size; ++n) {
		func.push_back(exp(2 * grid[n]));
		deriv.push_back(2 * exp(2 * grid[n]));
	}

	std::vector<double> res = mat_vec_product(3, deriv_operator, func);

	for (int n = 0; n != size; ++n) {
		std::cout 
			<< grid[n] << "\t" 
			<< func[n] << "\t"
			<< deriv[n] << "\t"
			<< res[n] << std::endl;
	}

	int int_cin;
	std::cin >> int_cin;

	return 0;

}