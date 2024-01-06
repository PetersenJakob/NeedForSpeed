#include "band_diagonal_matrix.h"


int main() {

	const int order = 6;

	TriDiagonal tri(order);
	PentaDiagonal penta(order);
	BandDiagonal band(order, 2, 1);

	print_matrix(tri);

	print_matrix(penta);

	print_matrix(band);

	return 0;
}