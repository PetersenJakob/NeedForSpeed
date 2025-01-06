#pragma once

#include <vector>


// Band-diagonal matrix stored in compact form.
template<typename T>
class BandDiagonalTemplate {

protected:

	// Matrix order: Number of elements along main diagonal.
	const std::size_t order_;
	// Lower bandwidth: Number of sub-diagonals.
	const std::size_t lower_bandwidth_;
	// Upper bandwidth: Number of super-diagonals.
	const std::size_t upper_bandwidth_;
	// Bandwidth: Maximum of lower_bandwidth_ and upper_bandwidth_.
	const std::size_t bandwidth_;
	// Number of diagonals.
	const std::size_t n_diagonals_;
	// Number of boundary rows (at each boundary).
	const std::size_t n_boundary_rows_;
	// Number of non-zero elements along each boundary row.
	const std::size_t n_boundary_elements_;

public:

	// Band-diagonal matrix in compact form.
	// TODO: Store as row-major or column-major? See Scott Meyers youtube video.
	// TODO: Should this be private/protected?
	std::vector<std::vector<T>> matrix;

	// Boundary rows of band-diagonal matrix.
	std::vector<std::vector<T>> boundary_rows;

	// TODO: Default constructor removed since useless...

	BandDiagonalTemplate(
		const std::size_t _order,
		const std::size_t _lower_bandwidth,
		const std::size_t _upper_bandwidth,
		const std::size_t _n_boundary_rows,
		const std::size_t _n_boundary_elements);

	BandDiagonalTemplate(const BandDiagonalTemplate& mat);

	// TODO: Why is this inherited by derived classes?
	bool operator==(const BandDiagonalTemplate& m);

	int order() const {
		return order_;
	}

	int lower_bandwidth() const {
		return lower_bandwidth_;
	}

	int upper_bandwidth() const {
		return upper_bandwidth_;
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
#if false
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
#endif
};


// Tri-diagonal matrix stored in compact form.
// TODO: Assumed, in other places, to have one boundary row!
// TODO: Remove _n_boundary_rows from parameter list!
template<typename T>
class TriDiagonalTemplate : public BandDiagonalTemplate<T> {

public:

	TriDiagonalTemplate() {}

	TriDiagonalTemplate(
		const std::size_t _order,
		const std::size_t _n_boundary_rows = 1,
		const std::size_t _n_boundary_elements = 3) :
		BandDiagonalTemplate<T>(_order, 1, 1, _n_boundary_rows, _n_boundary_elements) {}

	TriDiagonalTemplate(const TriDiagonalTemplate& mat) : BandDiagonalTemplate<T>(mat) {};
#if false
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
#endif
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
