#pragma once

#include <vector>

/**
 * @brief 
*/
class Bond {

private:

	std::vector<double> event_grid;
	int maturity_idx;

public:

	Bond();

	double maturity() {

		return event_grid[maturity_idx];

	}

	/**
	 * @brief 
	 * @param spot 
	 * @return 
	*/
	double payoff(const double spot);

	double price(
		const double spot,
		const int event_idx);

	/**
	 * @brief 
	 * @param spot 
	 * @param event_idx 
	 * @return 
	*/
	double delta(
		const double spot,
		const int event_idx);

	double gamma(
		const double spot,
		const int event_idx);

	double theta(
		const double spot,
		const int event_idx);

};
