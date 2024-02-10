#pragma once

#include <vector>

#include "band_diagonal_matrix.h"


void tri_solver(TriDiagonal& matrix, std::vector<double>& column);


void penta_solver(PentaDiagonal& matrix, std::vector<double>& column);


void tridiagonal_matrix_solver(
	std::vector<double>& sub,
	std::vector<double>& main, 
	std::vector<double>& super, 
	std::vector<double>& column,
	std::vector<double>& vec_tmp);


void pentadiagonal_matrix_solver(
	std::vector<double>& sub_2,
	std::vector<double>& sub_1,
	std::vector<double>& main,
	std::vector<double>& super_1,
	std::vector<double>& super_2,
	std::vector<double>& column,
	std::vector<double>& sub_tmp,
	std::vector<double>& main_tmp,
	std::vector<double>& super_tmp,
	std::vector<double>& vec_tmp);
