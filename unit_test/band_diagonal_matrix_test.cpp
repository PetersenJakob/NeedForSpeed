#include "pch.h"


TEST(Operator, X) {

	const int n_points = 21;

	TriDiagonal iden = identity::tri(n_points);

	print_matrix(iden);

	TriDiagonal mat1 = d1dx1::equidistant::c2b1(n_points, 0.1);

	print_matrix(mat1);

	print_matrix(13 * iden + mat1);

	TriDiagonal mat2 = d1dx1::equidistant::c2b1(n_points, 0.1);

	EXPECT_TRUE(mat1 == mat2);

	mat2 = 2.0 * mat1;

	EXPECT_TRUE(mat1 == 0.5 * mat2);

	EXPECT_NEAR(2.0, 2.0, 0.01);

}