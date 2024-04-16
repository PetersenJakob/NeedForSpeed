#pragma once


namespace vasicek {

	/**
	 * @brief Calculate A-function.
	 * 
	 * See Andersen & Piterbarg (2010), Proposition 10.1.4.
	 * 
	 * @param time_1: Initial time.
	 * @param time_2: Time of maturity.
	 * @param kappa: Speed of mean reversion.
	 * @param theta: Mean reversion level.
	 * @param sigma: Volatility.
	 * 
	 * @return A-function.
	*/
	double a_func(
		const double time_1,
		const double time_2,
		const double kappa,
		const double theta,
		const double sigma);

	/**
	 * @brief Calculate 1st order derivative of A-function wrt time_1.
	 * 
	 * See Andersen & Piterbarg (2010), Proposition 10.1.4.
	 * 
	 * @param time_1: Initial time.
	 * @param time_2: Time of maturity.
	 * @param kappa: Speed of mean reversion.
	 * @param theta: Mean reversion level.
	 * @param sigma: Volatility.
	 *
	 * @return Time derivative of A-function.
	*/
	double dadt_func(
		const double time_1,
		const double time_2,
		const double kappa,
		const double theta,
		const double sigma);


	/**
	 * @brief Calculate B-function. 
	 * 
	 * See Andersen & Piterbarg (2010), Proposition 10.1.4.
	 * 
	 * @param time_1: Initial time.
	 * @param time_2: Time of maturity.
	 * @param kappa: Speed of mean reversion.
	 * 
	 * @return B-function.
	*/
	double b_func(
		const double time_1,
		const double time_2,
		const double kappa);

	/**
	 * @brief Calculate 1st order derivative of B-function wrt time_1.
	 * 
	 * See Andersen & Piterbarg (2010), Proposition 10.1.4.
	 * 
	 * @param time_1: Initial time.
	 * @param time_2: Time of maturity.
	 * @param kappa: Speed of mean reversion.
	 * 
	 * @return Time derivative of B-function.
	*/
	double dbdt_func(
		const double time_1,
		const double time_2,
		const double kappa);

}
