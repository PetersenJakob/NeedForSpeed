#pragma once

#include <vector>

#include "band_diagonal_matrix.h"


namespace solver {

	// Band-diagonal matrix equation solver.
	void band(
		BandDiagonal& matrix,
		std::vector<double>& column);

	// Tri-diagonal matrix equation solver.
	void tri(
		BandDiagonal& matrix,
		std::vector<double>& column);

	// Penta-diagonal matrix equation solver.
	void penta(
		BandDiagonal& matrix,
		std::vector<double>& column);

}


void tridiagonal_matrix_solver(
	const std::vector<double>& sub,
	const std::vector<double>& main,
	const std::vector<double>& super,
	std::vector<double>& column,
	std::vector<double>& vec_tmp);


void pentadiagonal_matrix_solver(
	const std::vector<double>& sub_2,
	const std::vector<double>& sub_1,
	const std::vector<double>& main,
	const std::vector<double>& super_1,
	const std::vector<double>& super_2,
	std::vector<double>& column,
	std::vector<double>& sub_tmp,
	std::vector<double>& main_tmp,
	std::vector<double>& super_tmp,
	std::vector<double>& vec_tmp);
