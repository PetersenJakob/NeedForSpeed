#define _USE_MATH_DEFINES
#include <cmath>
#include <stdexcept>
#include <vector>

#include "heat_equation.h"


// The 3-dimensional heat equation reads
//		dV/dt = k * (d^2V/dx^2 + d^2V/dy^2 + d^2V/dz^2),
// where k > 0 is the thermal diffusivity. Using separation-of-varibles
//		V(t,x,y,z) = T(t) * X(x) * Y(y) * Z(z),
// the heat equation becomes
//		T' / (k * T) = X'' / X + Y'' / Y + Z'' / Z = -(a_x^2 + a_y^2 + a_z^2).
// 
// Assuming x in [0, L_x], t > 0, and defining
//		a_x(i) = i * pi / L_x, i \in N,
// a particular solution for the x-component, i.e., X'' = -a_x^2 * X, is given by
//		X_i(x) = sin(a_x(i) * x),
// where Dirichlet boundary conditions are used
//		X(0) = X(L_x) = 0.
// 
// Similar particular solutions for the y- and z-components are given by
// 		Y_j(y) = sin(a_y(j) * y), j \in N,
// 		Z_k(z) = sin(a_z(k) * z), k \in N,
// where 
//		a_y(j) = j * pi / L_y,
//		a_z(k) = k * pi / L_z.
// 
// The corresponding particular solution for the t-component is given by
//		T_{i,j,k}(t) = exp(-c_{i,j,k} * t),
// where
//		c_{i,j,k} = k * (a_x(i)^2 + a_y(j)^2 + a_z(k)^2).
//
// The full solution reads
//		V(t,x,y,z) = 
//			\sum_{i,j,k} d_{i,j,k} * T_{i,j,k}(t) * X_i(x) * Y_j(y) * Z_k(z),
// where d_{i,j,k} = 
//		(2 / L_x) * (2 / L_y) * (2 / L_z) * 
//			int_0^L_x int_0^L_y int_0^L_z 
//				V(x,y,z) * sin(a_x(i) * x) * sin(a_y(j) * y) * sin(a_z(k) * z) 
//					dx dy dz,
//  


// Particular solution for a spatial component.
std::vector<double> heat_eq::solution_space(
	const std::vector<double>& grid,
	const int order) {

	std::vector<double> solution(grid.size(), 0.0);

	const double grid_length = grid.back() - grid.front();

	for (int i = 0; i != grid.size(); ++i) {
		solution[i] = sin(order * M_PI * grid[i] / grid_length);
	}

	return solution;

}


// Particular solution for time component (for a given spatial component).
double heat_eq::solution_time(
	const double time,
	const std::vector<double>& grid,
	const int order,
	const double diffusivity) {

	const double grid_length = grid.back() - grid.front();

	const double factor = pow(order * M_PI / grid_length, 2);

	return exp(-factor * diffusivity * time);

}


// Full solution.
std::vector<double> heat_eq::solution_full(
	const double time,
	const std::vector<std::vector<double>>& grid,
	const std::vector<std::vector<int>>& order,
	const std::vector<std::vector<double>>& prefactor,
	const double diffusivity) {

	int n_points = 1;
	for (int i = 0; i != grid.size(); ++i) {
		n_points *= (int)grid[i].size();
	}
	std::vector<double> solution(n_points, 0.0);

	std::vector<std::vector<double>> particular_solutions;
	for (int i = 0; i != grid.size(); ++i) {
		std::vector<double> inner(grid[i].size(), 0.0);
		particular_solutions.push_back(inner);
	}

	// Loop over sets of "solution order".
	for (int i = 0; i != order.size(); ++i) {

		for (int j = 0; j != grid.size(); ++j) {
			// Particlar solution in time.
			const double solution_time_tmp =
				solution_time(time, grid[j], order[i][j], diffusivity);
			// Particular solution in space.
			particular_solutions[j] = solution_space(grid[j], order[i][j]);
			// Full particular solution.
			for (int k = 0; k != grid[j].size(); ++k) {
				particular_solutions[j][k] *= prefactor[i][j] * solution_time_tmp;
			}
		}

		int index = 0;
		if (grid.size() == 1) {
			for (int x1 = 0; x1 != grid[0].size(); ++x1) {
				solution[index] +=
					particular_solutions[0][x1];
				++index;
			}
		}
		else if (grid.size() == 2) {
			for (int x1 = 0; x1 != grid[0].size(); ++x1) {
				for (int x2 = 0; x2 != grid[1].size(); ++x2) {
					solution[index] +=
						particular_solutions[0][x1] *
						particular_solutions[1][x2];
					++index;
				}
			}
		}
		else if (grid.size() == 3) {
			for (int x1 = 0; x1 != grid[0].size(); ++x1) {
				for (int x2 = 0; x2 != grid[1].size(); ++x2) {
					for (int x3 = 0; x3 != grid[2].size(); ++x3) {
						solution[index] +=
							particular_solutions[0][x1] *
							particular_solutions[1][x2] *
							particular_solutions[2][x3];
						++index;
					}
				}
			}
		}
		else if (grid.size() == 4) {
			for (int x1 = 0; x1 != grid[0].size(); ++x1) {
				for (int x2 = 0; x2 != grid[1].size(); ++x2) {
					for (int x3 = 0; x3 != grid[2].size(); ++x3) {
						for (int x4 = 0; x4 != grid[3].size(); ++x4) {
							solution[index] +=
								particular_solutions[0][x1] *
								particular_solutions[1][x2] *
								particular_solutions[2][x3] *
								particular_solutions[3][x4];
							++index;
						}
					}
				}
			}
		}
		else if (grid.size() == 5) {
			for (int x1 = 0; x1 != grid[0].size(); ++x1) {
				for (int x2 = 0; x2 != grid[1].size(); ++x2) {
					for (int x3 = 0; x3 != grid[2].size(); ++x3) {
						for (int x4 = 0; x4 != grid[3].size(); ++x4) {
							for (int x5 = 0; x5 != grid[4].size(); ++x5) {
								solution[index] +=
									particular_solutions[0][x1] *
									particular_solutions[1][x2] *
									particular_solutions[2][x3] *
									particular_solutions[3][x4] *
									particular_solutions[4][x5];
								++index;
							}
						}
					}
				}
			}
		}
		else {
			throw std::invalid_argument("Dimension should be less than 6.");
		}

	}

	return solution;

}


std::function<std::vector<double>
	(const double, const std::vector<std::vector<double>>&)>
	heat_eq::solution_func(
		const std::vector<std::vector<int>>& order,
		const std::vector<std::vector<double>>& prefactor,
		const double diffusivity) {

	return [order, prefactor, diffusivity](
		const double time_,
		const std::vector<std::vector<double>>& grid_) {
			return heat_eq::solution_full(time_, grid_, order, prefactor, diffusivity);
		};

}
