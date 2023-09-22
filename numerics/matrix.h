#pragma once

#include <vector>

class Matrix {

public:
	Matrix();
	Matrix(int, int);

private:
	int n_rows;
	int n_cols;

	std::vector<std::vector<double>> data;

};
