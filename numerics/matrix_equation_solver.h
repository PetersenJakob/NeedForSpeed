#pragma once

#include <vector>

#include "band_diagonal_matrix.h"


namespace solver {

	// Tri-diagonal matrix equation solver.
	void tri(
		TriDiagonal& matrix, 
		std::vector<double>& column);

	// Penta-diagonal matrix equation solver.
	void penta(
		PentaDiagonal& matrix, 
		std::vector<double>& column);

	// Tri-diagonal matrix equation solver.
	void tri_test(
		BandDiagonal& matrix,
		std::vector<double>& column);

	// Penta-diagonal matrix equation solver.
	void penta_test(
		BandDiagonal& matrix,
		std::vector<double>& column);

}


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
