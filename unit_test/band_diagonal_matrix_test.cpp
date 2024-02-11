#include "pch.h"


TEST(TriDiagonal, Operators) {

	const bool print_data = false;

	const int n_points = 21;

	TriDiagonal iden = identity::tri(n_points);

	TriDiagonal tri1 = d1dx1::uniform::c2b1(n_points, 0.1);

	TriDiagonal tri2 = d1dx1::uniform::c2b1(n_points, 0.05);

	if (print_data) {

		print_matrix(iden);

		print_matrix(tri1);

		print_matrix(tri2);

		print_matrix(3 * iden + tri1);

	}

	EXPECT_TRUE(iden == iden);

	EXPECT_TRUE(2.0 * tri1 == tri1 * 2);

	EXPECT_TRUE(2 * iden + tri1 * 2.0 == 2 * (iden + tri1));

	EXPECT_TRUE(2.0 * tri1 == tri2);

	tri2 *= 0.5;

	EXPECT_TRUE(tri1 == tri2);

}


TEST(PentaDiagonal, Operators) {

	const bool print_data = false;

	const int n_points = 21;

	PentaDiagonal iden = identity::penta(n_points);

	std::vector<double> grid1 = grid::uniform(0.0, 1.0, n_points);
	std::vector<double> grid2 = grid::uniform(0.0, 0.5, n_points);

	PentaDiagonal penta1 = d1dx1::nonuniform::c4b2(n_points, grid1);

	PentaDiagonal penta2 = d1dx1::nonuniform::c4b2(n_points, grid2);

	if (print_data) {

		print_matrix(iden);

		print_matrix(penta1);

		print_matrix(penta2);

		print_matrix(3 * iden + penta1);

	}

	EXPECT_TRUE(iden == iden);

	EXPECT_TRUE(2.0 * penta1 == penta1 * 2);

	EXPECT_TRUE(2 * iden + penta1 * 2.0 == 2 * (iden + penta1));

	EXPECT_TRUE(2.0 * penta1 == penta2);

	penta2 *= 0.5;

	EXPECT_TRUE(penta1 == penta2);

}
