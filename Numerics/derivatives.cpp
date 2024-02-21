#include <stdexcept>
#include <vector>

#include "derivatives.h"

std::vector<element> d1dx1::central_2o_3p(int size, double dx) {
	
	// Total number of non-zero elements.
	int n_elements = size + 2 * (size - 1);

	// Container of non-zero elements.
	std::vector<element> container;
	container.reserve(n_elements);

	// Elements along sub-diagonal.
	for (int n = 1; n != size; ++n)
		container.push_back({ n, n - 1, 1.0 / 2.0 / dx });

	// Elements along main diagonal.
	for (int n = 0; n != size; ++n)
		container.push_back({ n, n, 0.0 });

	// Elements along super-diagonal.
	for (int n = 1; n != size; ++n)
		container.push_back({ n - 1, n, -1.0 / 2.0 / dx });

	return container;

}

std::vector<element> d1dx1::central_4o_5p(int size, double dx) {

	// Total number of non-zero elements.
	int n_elements = size + 2 * (size - 1) + 2 * (size - 2);

	// Container of non-zero elements.
	std::vector<element> container;
	container.reserve(n_elements);

	// Elements along second sub-diagonal.
	for (int n = 2; n != size; ++n)
		container.push_back({ n, n - 2, 1.0 / 2.0 / dx });

	// Elements along first sub-diagonal.
	for (int n = 1; n != size; ++n)
		container.push_back({ n, n - 1, 1.0 / 2.0 / dx });

	// Elements along main diagonal.
	for (int n = 0; n != size; ++n)
		container.push_back({ n, n, 0.0 });

	// Elements along first super-diagonal.
	for (int n = 1; n != size; ++n)
		container.push_back({ n - 1, n, -1.0 / 2.0 / dx });

	// Elements along second super-diagonal.
	for (int n = 2; n != size; ++n)
		container.push_back({ n - 2, n, -1.0 / 2.0 / dx });

	return container;

}

std::vector<double> mat_vec_product(const int stensil, const std::vector<element>& mat, const std::vector<double>& vec) {

	if (stensil != 3 and stensil != 5)
		throw std::invalid_argument("Stensil should be 3 or 5");

	int counter = 0;

	std::vector<double> res(vec.size(), 0.0);

	if (stensil == 5) {
		// Elements along second sub-diagonal.
		for (int n = 2; n != vec.size(); ++n) {
			res[n] += mat[counter].value * vec[n - 2];
			++counter;
		}
	}

	// Elements along first sub-diagonal.
	for (int n = 1; n != vec.size(); ++n) {
		res[n] += mat[counter].value * vec[n - 1];
		++counter;
	}

	// Elements along main diagonal.
	for (int n = 0; n != vec.size(); ++n) {
		res[n] += mat[counter].value * vec[n];
		++counter;
	}

	// Elements along first super-diagonal.
	for (int n = 1; n != vec.size(); ++n) {
		res[n - 1] += mat[counter].value * vec[n];
		++counter;
	}

	if (stensil == 5) {
		// Elements along second super-diagonal.
		for (int n = 2; n != vec.size(); ++n) {
			res[n - 2] += mat[counter].value * vec[n];
			++counter;
		}
	}
	
	return res;

}
