#pragma once


namespace vasicek {

	/**
	 * @brief
	 * @param time_1
	 * @param time_2
	 * @param kappa
	 * @return
	*/
	double b_func(
		const double time_1,
		const double time_2,
		const double kappa);


	double dbdt_func(
		const double time_1,
		const double time_2,
		const double kappa);


	double a_func(
		const double time_1,
		const double time_2,
		const double kappa,
		const double theta,
		const double sigma);


	double dadt_func(
		const double time_1,
		const double time_2,
		const double kappa,
		const double theta,
		const double sigma);

}
