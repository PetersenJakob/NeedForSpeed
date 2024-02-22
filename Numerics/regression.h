#pragma once

#include <vector>


namespace regression {

	// Simple linear regression.
	std::vector<double> slr(
		const std::vector<double>& x, 
		const std::vector<double>& y);

}