#pragma once

#include <Eigen/Dense>
#include <vector>

std::vector<double> grid_equidistant(const double x_min, const double x_max, const int n_points);

Eigen::VectorXd grid_eigen(const double x_min, const double x_max, const int n_points);
