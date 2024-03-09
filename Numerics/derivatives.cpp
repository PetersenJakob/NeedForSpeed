#include <vector>

#include "band_diagonal_matrix.h"
#include "coefficients.h"
#include "derivatives.h"
#include "utility.h"


// First order derivative operator. 
// Central difference; 2nd order accuracy. Boundary; 1st order accuracy.
TriDiagonal d1dx1::uniform::c2b1(const std::vector<double>& grid) {

	// TODO: If you choose two boundary rows with f1+f1 and b1+b1, what will the L2 function norm be?

	const int order = (int)grid.size();
	const double dx = grid[1] - grid[0];

	TriDiagonal matrix = setup<TriDiagonal>(order, coef_x1::uniform::c2(dx), 1, 2);
	boundary<TriDiagonal>(0, coef_x1::uniform::f1(dx), matrix);
	boundary<TriDiagonal>(1, coef_x1::uniform::b1(dx), matrix);

	return matrix;

}


// First order derivative operator. 
// Central difference; 2nd order accuracy. Boundary; 2nd order accuracy.
TriDiagonal d1dx1::uniform::c2b2(const std::vector<double>& grid) {

	const int order = (int)grid.size();
	const double dx = grid[1] - grid[0];

	TriDiagonal matrix = setup<TriDiagonal>(order, coef_x1::uniform::c2(dx), 1, 3);
	boundary<TriDiagonal>(0, coef_x1::uniform::f2(dx), matrix);
	boundary<TriDiagonal>(1, coef_x1::uniform::b2(dx), matrix);

	return matrix;

}


// First order derivative operator. 
// Central difference; 4th order accuracy. Boundary; 2nd order accuracy.
PentaDiagonal d1dx1::uniform::c4b2(const std::vector<double>& grid) {

	const int order = (int)grid.size();
	const double dx = grid[1] - grid[0];

	PentaDiagonal matrix = setup<PentaDiagonal>(order, coef_x1::uniform::c4(dx), 2, 3);
	boundary<PentaDiagonal>(0, coef_x1::uniform::f2(dx), matrix);
	boundary<PentaDiagonal>(1, coef_x1::uniform::f2(dx), matrix);
	boundary<PentaDiagonal>(2, coef_x1::uniform::b2(dx), matrix);
	boundary<PentaDiagonal>(3, coef_x1::uniform::b2(dx), matrix);

	return matrix;

}


// First order derivative operator. 
// Central difference; 4th order accuracy. Boundary; 4th order accuracy.
PentaDiagonal d1dx1::uniform::c4b4(const std::vector<double>& grid) {

	const int order = (int)grid.size();
	const double dx = grid[1] - grid[0];

	PentaDiagonal matrix = setup<PentaDiagonal>(order, coef_x1::uniform::c4(dx), 2, 5);
	boundary<PentaDiagonal>(0, coef_x1::uniform::f4(dx), matrix);
	boundary<PentaDiagonal>(1, coef_x1::uniform::f4(dx), matrix);
	boundary<PentaDiagonal>(2, coef_x1::uniform::b4(dx), matrix);
	boundary<PentaDiagonal>(3, coef_x1::uniform::b4(dx), matrix);

	return matrix;

}


// Second order derivative operator.
// Central difference; 2nd order accuracy. Boundary; d2dx2 = 0.
TriDiagonal d2dx2::uniform::c2b0(const std::vector<double>& grid) {

	const int order = (int)grid.size();
	const double dx = grid[1] - grid[0];

	TriDiagonal matrix = setup<TriDiagonal>(order, coef_x2::uniform::c2(dx), 1, 2);
	boundary<TriDiagonal>(0, { 0.0, 0.0 }, matrix);
	boundary<TriDiagonal>(1, { 0.0, 0.0 }, matrix);

	return matrix;

}


// Second order derivative operator.
// Central difference; 2nd order accuracy. Boundary; 1st order accuracy.
TriDiagonal d2dx2::uniform::c2b1(const std::vector<double>& grid) {

	const int order = (int)grid.size();
	const double dx = grid[1] - grid[0];

	TriDiagonal matrix = setup<TriDiagonal>(order, coef_x2::uniform::c2(dx), 1, 3);
	boundary<TriDiagonal>(0, coef_x2::uniform::f1(dx), matrix);
	boundary<TriDiagonal>(1, coef_x2::uniform::b1(dx), matrix);

	return matrix;

}


// Second order derivative operator.
// Central difference; 2nd order accuracy. Boundary; 2nd order accuracy.
TriDiagonal d2dx2::uniform::c2b2(const std::vector<double>& grid) {

	const int order = (int)grid.size();
	const double dx = grid[1] - grid[0];

	TriDiagonal matrix = setup<TriDiagonal>(order, coef_x2::uniform::c2(dx), 1, 4);
	boundary<TriDiagonal>(0, coef_x2::uniform::f2(dx), matrix);
	boundary<TriDiagonal>(1, coef_x2::uniform::b2(dx), matrix);

	return matrix;

}


// Second order derivative operator.
// Central difference; 4th order accuracy. Boundary; 2nd row c2, 1st row d2dx2 = 0.
PentaDiagonal d2dx2::uniform::c4b0(const std::vector<double>& grid) {

	const int order = (int)grid.size();
	const double dx = grid[1] - grid[0];

	PentaDiagonal matrix = setup<PentaDiagonal>(order, coef_x2::uniform::c4(dx), 2, 3);

	boundary<PentaDiagonal>(0, { 0.0, 0.0, 0.0 }, matrix);


	// Central difference at 2nd row.
	matrix.boundary_rows[1][0] = coef_x2::uniform::c2(dx)[0];
	matrix.boundary_rows[1][1] = coef_x2::uniform::c2(dx)[1];
	matrix.boundary_rows[1][2] = coef_x2::uniform::c2(dx)[2];

	matrix.boundary_rows[2][1] = coef_x2::uniform::c2(dx)[0];
	matrix.boundary_rows[2][2] = coef_x2::uniform::c2(dx)[1];
	matrix.boundary_rows[2][3] = coef_x2::uniform::c2(dx)[2];


	boundary<PentaDiagonal>(3, { 0.0, 0.0, 0.0 }, matrix);

	return matrix;

}


// Second order derivative operator.
// Central difference; 4th order accuracy. Boundary; TODO 1st order accuracy.
PentaDiagonal d2dx2::uniform::c4b4(const std::vector<double>& grid) {

	const int order = (int)grid.size();
	const double dx = grid[1] - grid[0];

#if false
	PentaDiagonal matrix = setup<PentaDiagonal>(order, coef_x2::uniform::c4(dx), 2, 3);
	boundary<PentaDiagonal>(0, coef_x2::uniform::f1(dx), matrix);
	boundary<PentaDiagonal>(1, coef_x2::uniform::f1(dx), matrix);
	boundary<PentaDiagonal>(2, coef_x2::uniform::b1(dx), matrix);
	boundary<PentaDiagonal>(3, coef_x2::uniform::b1(dx), matrix);
#endif

	PentaDiagonal matrix = setup<PentaDiagonal>(order, coef_x2::uniform::c4(dx), 2, 6);
	boundary<PentaDiagonal>(0, coef_x2::uniform::f4(dx), matrix);
	boundary<PentaDiagonal>(1, coef_x2::uniform::f4(dx), matrix);
	boundary<PentaDiagonal>(2, coef_x2::uniform::b4(dx), matrix);
	boundary<PentaDiagonal>(3, coef_x2::uniform::b4(dx), matrix);

	return matrix;

}


// First order derivative operator.
// Central difference; 2nd order accuracy. Boundary; 1st order accuracy.
//TriDiagonal d1dx1::nonuniform::c2b1(const int order, const std::vector<double> grid) {
TriDiagonal d1dx1::nonuniform::c2b1(const std::vector<double>& grid) {

	const int order = (int)grid.size();

	TriDiagonal matrix = setup<TriDiagonal>(order, grid, coef_x1::nonuniform::c2, 1, 2);

	std::vector<double> dx_vec_first = { 0.0, 0.0, grid[1] - grid[0], 0.0 };
	boundary<TriDiagonal>(0, coef_x1::nonuniform::f1(dx_vec_first), matrix);
	// TODO: Should dx_last be reversed?
	std::vector<double> dx_vec_last = { 0.0, 0.0, grid[order - 1] - grid[order - 2], 0.0 };
	boundary<TriDiagonal>(1, coef_x1::nonuniform::b1(dx_vec_last), matrix);

	return matrix;

}


// First order derivative operator.
// Central difference; 2nd order accuracy. Boundary; 2nd order accuracy.
//TriDiagonal d1dx1::nonuniform::c2b2(const int order, const std::vector<double> grid) {
TriDiagonal d1dx1::nonuniform::c2b2(const std::vector<double>& grid) {

	const int order = (int)grid.size();

	TriDiagonal matrix = setup<TriDiagonal>(order, grid, coef_x1::nonuniform::c2, 1, 3);

	std::vector<double> dx_vec_first = { 0.0, 0.0, grid[1] - grid[0], grid[2] - grid[1] };
	boundary<TriDiagonal>(0, coef_x1::nonuniform::f2(dx_vec_first), matrix);

	// TODO: Should dx_last be reversed? Yes!
	std::vector<double> dx_vec_last = { 0.0, 0.0, grid[order - 1] - grid[order - 2], grid[order - 2] - grid[order - 3] };
	boundary<TriDiagonal>(1, coef_x1::nonuniform::b2(dx_vec_last), matrix);

	return matrix;

}


// First order derivative operator.
// Central difference; 4th order accuracy. Boundary; 2nd order accuracy.
PentaDiagonal d1dx1::nonuniform::c4b2(const std::vector<double>& grid) {

	// TODO: Be able to choose c2 as second boundary row, and f1 as first boundary row.

	const int order = (int)grid.size();

	PentaDiagonal matrix = setup<PentaDiagonal>(order, grid, coef_x1::nonuniform::c4, 2, 3);

	std::vector<double> dx_vec_1 = { 0.0, 0.0, grid[1] - grid[0], grid[2] - grid[1] };
	boundary<PentaDiagonal>(0, coef_x1::nonuniform::f2(dx_vec_1), matrix);

	std::vector<double> dx_vec_2 = { 0.0, 0.0, grid[2] - grid[1], grid[3] - grid[2] };
	boundary<PentaDiagonal>(1, coef_x1::nonuniform::f2(dx_vec_2), matrix);

	std::vector<double> dx_vec_3 = { 0.0, 0.0, grid[order - 2] - grid[order - 3], grid[order - 3] - grid[order - 4] };
	boundary<PentaDiagonal>(2, coef_x1::nonuniform::b2(dx_vec_3), matrix);

	std::vector<double> dx_vec_4 = { 0.0, 0.0, grid[order - 1] - grid[order - 2], grid[order - 2] - grid[order - 3] };
	boundary<PentaDiagonal>(3, coef_x1::nonuniform::b2(dx_vec_4), matrix);

	return matrix;

}


// Second order derivative operator.
// Central difference; 2nd order accuracy. Boundary; d2dx2 = 0.
TriDiagonal d2dx2::nonuniform::c2b0(const std::vector<double>& grid) {

	const int order = (int)grid.size();

	TriDiagonal matrix = setup<TriDiagonal>(order, grid, coef_x2::nonuniform::c2, 1, 2);

	boundary<TriDiagonal>(0, { 0.0, 0.0 }, matrix);
	boundary<TriDiagonal>(1, { 0.0, 0.0 }, matrix);

	return matrix;

}


// Second order derivative operator.
// Central difference; 2nd order accuracy. Boundary; 1st order accuracy.
TriDiagonal d2dx2::nonuniform::c2b1(const std::vector<double>& grid) {

	const int order = (int)grid.size();

	TriDiagonal matrix = setup<TriDiagonal>(order, grid, coef_x2::nonuniform::c2, 1, 3);

	std::vector<double> dx_vec_first = { 0.0, 0.0, grid[1] - grid[0], grid[2] - grid[1] };
	boundary<TriDiagonal>(0, coef_x2::nonuniform::f1(dx_vec_first), matrix);
	// TODO: Should dx_last be reversed?
	std::vector<double> dx_vec_last = { 0.0, 0.0, grid[order - 1] - grid[order - 2], grid[order - 2] - grid[order - 3] };
	boundary<TriDiagonal>(1, coef_x2::nonuniform::b1(dx_vec_last), matrix);

	return matrix;

}


// Second order derivative operator.
// Central difference; 4th order accuracy. Boundary; 2nd row c2, 1st row d2dx2 = 0.
PentaDiagonal d2dx2::nonuniform::c4b0(const std::vector<double>& grid) {

	const int order = (int)grid.size();

	PentaDiagonal matrix = setup<PentaDiagonal>(order, grid, coef_x2::nonuniform::c4, 2, 3);

	std::vector<double> dx_vec_1 = { 0.0, 0.0, 0.0 };
	boundary<PentaDiagonal>(0, dx_vec_1, matrix);


	// Central difference at 2nd row.
	std::vector<double> dx_vec_2 = { 0.0, grid[1] - grid[0], grid[2] - grid[1], 0.0 };
	matrix.boundary_rows[1][0] = coef_x2::nonuniform::c2(dx_vec_2)[0];
	matrix.boundary_rows[1][1] = coef_x2::nonuniform::c2(dx_vec_2)[1];
	matrix.boundary_rows[1][2] = coef_x2::nonuniform::c2(dx_vec_2)[2];

	std::vector<double> dx_vec_3 = { 0.0, grid[order - 2] - grid[order - 3], grid[order - 1] - grid[order - 2], 0.0 };
	matrix.boundary_rows[2][1] = coef_x2::nonuniform::c2(dx_vec_3)[0];
	matrix.boundary_rows[2][2] = coef_x2::nonuniform::c2(dx_vec_3)[1];
	matrix.boundary_rows[2][3] = coef_x2::nonuniform::c2(dx_vec_3)[2];


	std::vector<double> dx_vec_4 = { 0.0, 0.0, 0.0 };
	boundary<PentaDiagonal>(3, dx_vec_4, matrix);

	return matrix;

}
