#pragma once

#include <stdexcept>
#include <vector>


// TODO: POSSIBLE LINING ERRORS SINCE ALL CLASS METHODS SHOULD BE DEFINED IN THIS HEADER!

// Band-diagonal matrix stored in compact form.
// TODO: General assumption. For non-boundary rows, all diagonal elementss are "present"!!!
// TODO: This assumption has been used in matrix_multiply_column function!
// TODO: order- >= n_boundary_rows_lower + n_boundary_rows_upper, included assert in constructor.
// TODO: Upper boundary row vectors are intepreted as being from right to left
template<typename Tnumber>
class BandDiagonalTemplate {

protected:

	// TODO: Make "dim" struct with all data members?

	// Matrix order: Number of elements along main diagonal.
	const std::size_t order_;
	// Lower bandwidth: Number of sub-diagonals.
	const std::size_t bandwidth_lower_;
	// Upper bandwidth: Number of super-diagonals.
	const std::size_t bandwidth_upper_;
	// Bandwidth: Maximum of bandwidth_lower_ and bandwidth_upper_.
	const std::size_t bandwidth_;
	// Number of diagonals.
	const std::size_t n_diagonals_;
	// Number of rows at lower boundary.
	const std::size_t n_boundary_rows_lower_;
	// Number of rows at upper boundary.
	const std::size_t n_boundary_rows_upper_;
	// Number of boundary rows.
	const std::size_t n_boundary_rows_;
	// Number of non-zero elements along each lower boundary row.
	const std::size_t n_boundary_elements_lower_;
	// Number of non-zero elements along each upper boundary row.
	const std::size_t n_boundary_elements_upper_;

public:

	// Band-diagonal matrix in compact form.
	// TODO: Store as row-major or column-major? See Scott Meyers youtube video.
	// TODO: Should this be private/protected?
	std::vector<std::vector<Tnumber>> matrix;
	// Lower boundary rows of band-diagonal matrix.
	std::vector<std::vector<Tnumber>> boundary_lower;
	// Upper boundary rows of band-diagonal matrix.
	std::vector<std::vector<Tnumber>> boundary_upper;

	// TODO: Default constructor removed since useless...

	BandDiagonalTemplate(
		const std::size_t _order,
		const std::size_t _bandwidth_lower,
		const std::size_t _bandwidth_upper,
		const std::size_t _n_boundary_rows_lower,
		const std::size_t _n_boundary_rows_upper,
		const std::size_t _n_boundary_elements_lower,
		const std::size_t _n_boundary_elements_upper);

	BandDiagonalTemplate(const BandDiagonalTemplate& mat);

	// TODO: Why is this inherited by derived classes?
	bool operator==(const BandDiagonalTemplate& rhs);

	BandDiagonalTemplate operator+(const BandDiagonalTemplate& rhs);

	BandDiagonalTemplate& operator+=(const BandDiagonalTemplate& rhs);

	// TODO: const reference? Tnumber object might be "big"?
	BandDiagonalTemplate operator+(const Tnumber rhs);

	// TODO: const reference? Tnumber object might be "big"?
	BandDiagonalTemplate& operator+=(const Tnumber rhs);

	BandDiagonalTemplate operator-(const BandDiagonalTemplate& rhs);

	BandDiagonalTemplate& operator-=(const BandDiagonalTemplate& rhs);

	// TODO: const reference? Tnumber object might be "big"?
	BandDiagonalTemplate operator-(const Tnumber rhs);

	// TODO: const reference? Tnumber object might be "big"?
	BandDiagonalTemplate& operator-=(const Tnumber rhs);

	// TODO: const reference? Tnumber object might be "big"?
	BandDiagonalTemplate operator*(const Tnumber rhs);

	// TODO: const reference? Tnumber object might be "big"?
	BandDiagonalTemplate& operator*=(const Tnumber rhs);

	std::vector<Tnumber> operator*(const std::vector<Tnumber>& rhs);

	int order() const {
		return order_;
	}

	int bandwidth_lower() const {
		return bandwidth_lower_;
	}

	int bandwidth_upper() const {
		return bandwidth_upper_;
	}

	int bandwidth() const {
		return bandwidth_;
	}

	int n_diagonals() const {
		return n_diagonals_;
	}

	int n_boundary_rows_lower() const {
		return n_boundary_rows_lower_;
	}

	int n_boundary_rows_upper() const {
		return n_boundary_rows_upper_;
	}

	int n_boundary_rows() const {
		return n_boundary_rows_;
	}

	int n_boundary_elements_lower() const {
		return n_boundary_elements_lower_;
	}

	int n_boundary_elements_upper() const {
		return n_boundary_elements_upper_;
	}

	BandDiagonalTemplate identity();

	BandDiagonalTemplate pre_vector(const std::vector<Tnumber>& vector);

#if false
	// ...
	template<typename Tnumber>
	void gauss_elimination(
		const std::size_t boundary_row_idx,
		const std::size_tboundary_element_idx,
		const std::size_t matrix_row_idx,
		std::vector<Tnumber>& column) {};

	// ...
	void overwrite_bounary_row(const std::size_t boundary_row_idx) {};

	// ...
	template<typename Tnumber>
	virtual void adjust_boundary(std::vector<Tnumber>& column) {};

	void print_matrixTemplate() {};

#endif
};


template<typename Tnumber>
BandDiagonalTemplate<Tnumber>::BandDiagonalTemplate(
	const std::size_t _order,
	const std::size_t _bandwidth_lower,
	const std::size_t _bandwidth_upper,
	const std::size_t _n_boundary_rows_lower,
	const std::size_t _n_boundary_rows_upper,
	const std::size_t _n_boundary_elements_lower,
	const std::size_t _n_boundary_elements_upper) {

	order_ = _order;
	bandwidth_lower_ = _bandwidth_lower;
	bandwidth_upper_ = _bandwidth_upper;
	bandwidth_ = std::max(bandwidth_lower_, bandwidth_upper_);
	n_diagonals_ = 1 + bandwidth_lower_ + bandwidth_upper_;
	n_boundary_rows_lower_ = _n_boundary_rows_lower;
	n_boundary_rows_upper_ = _n_boundary_rows_upper;
	n_boundary_rows_ = n_boundary_rows_lower_ + n_boundary_rows_upper_;
	n_boundary_elements_lower_ = _n_boundary_elements_lower;
	n_boundary_elements_upper_ = _n_boundary_elements_upper;

	if (order_ < n_boundary_rows_lower_ + n_boundary_rows_upper_) {
		throw std::invalid_argument("order < n_boundary_rows_lower + n_boundary_rows_upper!");
	}

	// Band-diagonal matrix in compact form.
	std::vector<Tnumber> diagonal(order_, Tnumber(0.0));
	std::vector<std::vector<Tnumber>> m(n_diagonals_, diagonal);
	matrix = std::move(m);

	// Lower boundary rows of band-diagonal matrix.
	std::vector<Tnumber> row_lower(n_boundary_elements_lower_, Tnumber(0.0));
	std::vector<std::vector<Tnumber>> b_lower(2 * n_boundary_rows_lower_, row_lower);
	boundary_lower = std::move(b_lower);

	// Upper boundary rows of band-diagonal matrix.
	std::vector<Tnumber> row_upper(n_boundary_elements_upper_, Tnumber(0.0));
	std::vector<std::vector<Tnumber>> b_upper(2 * n_boundary_rows_upper_, row_upper);
	boundary_upper = std::move(b_upper);

}


template<typename Tnumber>
BandDiagonalTemplate<Tnumber>::BandDiagonalTemplate(const BandDiagonalTemplate& mat) {

	order_ = mat.order_;
	bandwidth_lower_ = mat.bandwidth_lower_;
	bandwidth_upper_ = mat.bandwidth_upper_;
	bandwidth_ = mat.bandwidth_;
	n_diagonals_ = mat.n_diagonals_;
	n_boundary_rows_lower_ = mat.n_boundary_rows_lower_;
	n_boundary_rows_upper_ = mat.n_boundary_rows_upper_;
	n_boundary_rows_ = mat.n_boundary_rows_;
	n_boundary_elements_lower_ = mat.n_boundary_elements_lower_;
	n_boundary_elements_upper_ = mat.n_boundary_elements_upper_;
	matrix = mat.matrix;
	boundary_lower = mat.boundary_lower;
	boundary_upper = mat.boundary_upper;

}


template<typename Tnumber>
void matrix_add_matrixTemplate(
	const BandDiagonalTemplate<Tnumber>& matrix,
	BandDiagonalTemplate<Tnumber>& result) {

	// TODO: Check if the two matrices have same dimensions?

	// TODO: Maybe adjust when row-major/column-major structure has been chosen?
	for (std::size_t i = 0; i != result.n_diagonals(); ++i) {
		for (std::size_t j = 0; j != result.order(); ++j) {
			result.matrix[i][j] += matrix.matrix[i][j];
		}
	}

	for (std::size_t i = 0; i != result.n_boundary_rows_lower(); ++i) {
		for (std::size_t j = 0; j != result.n_boundary_elements_lower(); ++j) {
			result.boundary_lower[i][j] += matrix.boundary_lower[i][j];
		}
	}

	for (std::size_t i = 0; i != result.n_boundary_rows_upper(); ++i) {
		for (std::size_t j = 0; j != result.n_boundary_elements_upper(); ++j) {
			result.boundary_upper[i][j] += matrix.boundary_upper[i][j];
		}
	}

}


template<typename Tnumber>
void matrix_add_scalarTemplate(
	BandDiagonalTemplate<Tnumber>& matrix,
	const Tnumber scalar) {

	for (std::size_t i = 0; i != matrix.n_diagonals(); ++i) {
		for (std::size_t j = 0; j != matrix.order(); ++j) {
			matrix.matrix[i][j] += scalar;
		}
	}

	for (std::size_t i = 0; i != matrix.n_boundary_rows_lower(); ++i) {
		for (std::size_t j = 0; j != matrix.n_boundary_elements_lower(); ++j) {
			matrix.boundary_lower[i][j] = +scalar;
		}
	}

	for (std::size_t i = 0; i != matrix.n_boundary_rows_upper(); ++i) {
		for (std::size_t j = 0; j != matrix.n_boundary_elements_upper(); ++j) {
			matrix.boundary_upper[i][j] = +scalar;
		}
	}

}


template<typename Tnumber>
void matrix_multiply_scalarTemplate(
	BandDiagonalTemplate<Tnumber>& matrix,
	const Tnumber scalar) {

	// TODO: Maybe adjust when row-major/column-major structure has been chosen?
	for (std::size_t i = 0; i != matrix.n_diagonals(); ++i) {
		for (std::size_t j = 0; j != matrix.order(); ++j) {
			matrix.matrix[i][j] *= scalar;
		}
	}

	for (std::size_t i = 0; i != matrix.n_boundary_rows_lower(); ++i) {
		for (std::size_t j = 0; j != matrix.n_boundary_elements_lower(); ++j) {
			matrix.boundary_lower[i][j] *= scalar;
		}
	}

	for (std::size_t i = 0; i != matrix.n_boundary_rows_upper(); ++i) {
		for (std::size_t j = 0; j != matrix.n_boundary_elements_upper(); ++j) {
			matrix.boundary_upper[i][j] *= scalar;
		}
	}

}


template<typename Tnumber>
BandDiagonalTemplate<Tnumber> operator*(
	const Tnumber scalar,
	const BandDiagonalTemplate<Tnumber> rhs);


template<typename Tnumber>
void matrix_multiply_columnTemplate(
	const BandDiagonalTemplate<Tnumber>& matrix,
	const std::vector<Tnumber>& column,
	std::vector<Tnumber>& result) {

	std::size_t row_idx = 0;
	std::size_t column_idx = 0;

	// Contribution from lower boundary rows of matrix.
	for (int i = 0; i != matrix.n_boundary_rows_lower(); ++i) {
		for (int j = 0; j != matrix.n_boundary_elements_lower(); ++j) {
			result[i] += matrix.boundary_lower[i][j] * column[j];
		}
	}

	// Contribution from upper boundary rows of matrix.
	for (int i = 0; i != matrix.n_boundary_rows_upper(); ++i) {
		for (int j = 0; j != matrix.n_boundary_elements_upper(); ++j) {
			row_idx = (matrix.order() - 1) - i;
			column_idx = (matrix.order() - 1) - j;
			result[row_idx] += matrix.boundary_upper[i][j] * column[column_idx];
		}
	}

	// TODO: Only use interior rows from matrix. 
	// TODO: Elements of boundary rows should be represented by boundary_rows.

	// Contribution from "interior" rows of matrix.
	const std::size_t i_initial = matrix.n_boundary_rows_lower();
	const std::size_t i_final = matrix.order() - matrix.n_boundary_rows_upper();
	const std::size_t j_initial = 0;
	const std::size_t j_final = matrix.n_diagonals();

	for (std::size_t i = i_initial; i != i_final; ++i) {
		for (std::size_t j = j_initial; j != j_final; ++j) {
			result[i] += matrix.matrix[j][i]
				* column[(i - matrix.n_boundary_rows_lower()) + j];
		}
	}

}


template<typename Tnumber>
void prevector_multiply_matrixTemplate(
	const std::vector<Tnumber>& vector,
	BandDiagonalTemplate<Tnumber>& matrix);


// Tri-diagonal matrix stored in compact form.
// TODO: Assumed, in other places, to have one boundary row!
// TODO: Remove _n_boundary_rows from parameter list!
template<typename Tnumber>
class TriDiagonalTemplate : public BandDiagonalTemplate<Tnumber> {

public:

	TriDiagonalTemplate(
		const std::size_t _order,
		const std::size_t _n_boundary_elements_lower = 2,
		const std::size_t _n_boundary_elements_upper = 2) :
		BandDiagonalTemplate<Tnumber>(_order, 1, 1, 1, 1, 
			_n_boundary_elements_lower, _n_boundary_elements_upper) {}

	// TODO: Necessary?
	TriDiagonalTemplate(const TriDiagonalTemplate& mat) : 
		BandDiagonalTemplate<Tnumber>(mat) {};

};


// Penta-diagonal matrix stored in compact form.
// TODO: Assumed, in other places, to have two boundary row!
// TODO: Remove _n_boundary_rows from parameter list!
template<typename Tnumber>
class PentaDiagonalTemplate : public BandDiagonalTemplate<Tnumber> {

public:

	PentaDiagonalTemplate(
		const int _order,
		const int _n_boundary_elements_lower = 3,
		const int _n_boundary_elements_upper = 3) :
		BandDiagonalTemplate<Tnumber>(_order, 2, 2, 2, 2, 
			_n_boundary_elements_lower, _n_boundary_elements_upper) {}

	// TODO: Necessary?
	PentaDiagonalTemplate(const PentaDiagonalTemplate& mat) :
		BandDiagonalTemplate<Tnumber>(mat) {};

};


// ###############################################################################


// Band-diagonal matrix stored in compact form.
// Note: The lower and upper bandwidths are assumed to be identical.
class BandDiagonal {

protected:

	// Matrix order: Number of elements along main diagonal.
	int order_;
	// Bandwidth: Number of sub-diagonals or super-diagonals.
	int bandwidth_;
	// Number of diagonals.
	int n_diagonals_;
	// Number of boundary rows (at each boundary).
	int n_boundary_rows_;
	// Number of non-zero elements along each boundary row.
	int n_boundary_elements_;

public:

	// Band-diagonal matrix in compact form.
	std::vector<std::vector<double>> matrix;

	// Boundary rows of band-diagonal matrix.
	std::vector<std::vector<double>> boundary_rows;

	std::vector<std::vector<double>> boundary_rows_tmp;

	BandDiagonal();

	BandDiagonal(
		const int _order,
		const int _bandwidth,
		const int _n_boundary_rows,
		const int _n_boundary_elements);

	BandDiagonal(const BandDiagonal& mat);

	// Why is this inherited by derived classes?
	bool operator==(const BandDiagonal& m);

	int order() const {
		return order_;
	}

	int bandwidth() const {
		return bandwidth_;
	}

	int n_diagonals() const {
		return n_diagonals_;
	}

	int n_boundary_rows() const {
		return n_boundary_rows_;
	}

	int n_boundary_elements() const {
		return n_boundary_elements_;
	}

	// ...
	void gauss_elimination(
		const int boundary_row_idx,
		const int boundary_element_idx,
		const int matrix_row_idx,
		std::vector<double>& column);

	// ...
	void overwrite_bounary_row(const int boundary_row_idx);

	// ...
	virtual void adjust_boundary(std::vector<double>& column) {}

};


// Tri-diagonal matrix stored in compact form.
// TODO: Assumed, in other places, to have one boundary row!
// Remove _n_boundary_rows from parameter list!
class TriDiagonal : public BandDiagonal {

public:

	TriDiagonal() {}

	TriDiagonal(
		const int _order,
		const int _n_boundary_rows = 1,
		const int _n_boundary_elements = 3) : 
		BandDiagonal(_order, 1, _n_boundary_rows, _n_boundary_elements) {}

	TriDiagonal operator*(const double scalar);

	std::vector<double> operator*(const std::vector<double>& vector);

	TriDiagonal operator*=(const double scalar);

	TriDiagonal operator+(const TriDiagonal& rhs);

	TriDiagonal operator+=(const TriDiagonal& rhs);

	TriDiagonal operator-(const TriDiagonal& rhs);

	TriDiagonal operator-=(const TriDiagonal& rhs);

	TriDiagonal identity();

	TriDiagonal pre_vector(const std::vector<double>& vector);

	void adjust_boundary(std::vector<double>& column);

};


// Penta-diagonal matrix stored in compact form.
// TODO: Assumed, in other places, to have two boundary row!
// Remove _n_boundary_rows from parameter list!
class PentaDiagonal : public BandDiagonal {

public:

	PentaDiagonal() {}

	PentaDiagonal(
		const int _order,
		const int _n_boundary_rows = 2,
		const int _n_boundary_elements = 3) : 
		BandDiagonal(_order, 2, _n_boundary_rows, _n_boundary_elements) {}

	PentaDiagonal operator*(const double scalar);

	std::vector<double> operator*(const std::vector<double>& vector);

	PentaDiagonal operator*=(const double scalar);

	PentaDiagonal operator+(const PentaDiagonal& rhs);

	PentaDiagonal operator+=(const PentaDiagonal& rhs);

	PentaDiagonal operator-(const PentaDiagonal& rhs);

	PentaDiagonal operator-=(const PentaDiagonal& rhs);

	PentaDiagonal identity();

	PentaDiagonal pre_vector(const std::vector<double>& vector);

	void adjust_boundary(std::vector<double>& column);

};


TriDiagonal operator*(
	const double scalar, 
	TriDiagonal rhs);


PentaDiagonal operator*(
	const double scalar, 
	PentaDiagonal rhs);


template<class T>
void scalar_multiply_matrix(
	const double scalar, 
	T& matrix);


template<class T>
void matrix_multiply_vector(
	const T& matrix, 
	const std::vector<double>& vector, 
	std::vector<double>& result);


template<class T>
void matrix_add_matrix(
	const T& matrix1, 
	const T& matrix2, 
	T& result);


void print_matrix(BandDiagonal matrix);


template<class T>
void row_multiply_matrix(
	T& matrix,
	const std::vector<double>& vector);
