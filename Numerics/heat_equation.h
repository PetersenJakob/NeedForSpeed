#pragma once

#include <vector>


std::vector<double> solution_space(
	const std::vector<double>& grid,
	const int order);


double solution_time(
	const double time,
	const std::vector<double>& grid,
	const int order,
	const double diffusivity);


std::vector<double> solution_full(
	const double time,
	const std::vector<std::vector<double>>& grid,
	const std::vector<std::vector<int>>& order,
	const std::vector<std::vector<double>>& prefactor,
	const double diffusivity);
