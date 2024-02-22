#pragma once


#include <stdexcept>
#include <typeinfo>
#include <vector>

#include "band_diagonal_matrix.h"

template <class T>
void theta_propagator(
	const double dt,
	T& identity_operator,
	T& derivative_operator,
	std::vector<double>& column,
	const double theta_parameter = 0.5) {

	// Left-hand-side operator.
	T lhs = derivative_operator;
	lhs *= -theta_parameter * dt;
	lhs += identity_operator;

	// RHS operator.
	T rhs = derivative_operator;
	rhs *= (1.0 - theta_parameter) * dt;
	rhs += identity_operator;

	// Evaluation RHS.
	column = rhs * column;

	// Solve matrix equation.
	if (typeid(lhs) == typeid(TriDiagonal)) {
		solver::tri_test(lhs, column);
	}
	else if (typeid(lhs) == typeid(PentaDiagonal)) {
		solver::penta_test(lhs, column);
	}
	else {
		throw std::invalid_argument("BandDiagonal object unknown.");
	}

}
