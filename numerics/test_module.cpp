#include "band_diagonal_matrix.h"
#include "finite_difference.h"


int main() {

	const int order = 6;

	TriDiagonal tri(order);
	PentaDiagonal penta(order);
	BandDiagonal band(order, 2, 1);

	print_matrix(tri);

	const double dx = 0.1;
	print_matrix(d1dx1::c2(order, dx));

	print_matrix(penta);

//	print_matrix(band);

	return 0;
}