#pragma once

#include <vector>

#include "band_diagonal_matrix.h"


void tri_solver(TriDiagonal& matrix, std::vector<double>& column);


void tridiagonal_matrix_solver(
	std::vector<double>& sub,
	std::vector<double>& main, 
	std::vector<double>& super, 
	std::vector<double>& column);
