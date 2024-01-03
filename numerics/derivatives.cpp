#include <vector>

#include "derivatves.h"

std::vector<triple> central(int n_points) {
	
	std::vector<triple> derivative;

	// Sub-diagonal.

	// Main diagonal.
	for (int n = 0; n != n_points; n++) {
		derivative.push_back({ 1, 1, 0 });
	}

	// Super-diagonal.

	return derivative;

}