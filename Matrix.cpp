#include "Matrix.h"
#include<map>
#include<utility>
#include<initializer_list>
#include<stdexcept>
#include<string>
#include<sstream>
#include<cmath>
#include<functional>
#include<algorithm>
#include<iostream>

namespace matrix
{
	struct Matrix::_Matrix {
		//! An ordered pair of coordinates, for accessing values in \link #numbers\endlink.
		typedef std::pair<coord_t, coord_t> orderedPairType;
		//! The type of \link #numbers\endlink.
		typedef std::map<orderedPairType, double> numberMapType;

		//! Recursively calculates determinant of matrix.
		double determinant (const Matrix& matrix)
		{
			dimens_t sideLength = matrix.width();

			//Precondition: matrix is square
			if (sideLength == 1) {
				return matrix (1, 1);
			} else if (sideLength == 2) {
				return matrix (1, 1) * matrix (2, 2)
				       - matrix (1, 2) * matrix (2, 1);
			} else {
				std::vector<coord_t> coords (sideLength, 1);
				//Generate numbers [1,sideLength]
				coord_t n = 2;
				std::generate (coords.begin() + 1, coords.end(),
				               [&n] {return n++;});

				double determinant {};
				/* A sign function for compact sign changing.

				Since we'll always use row 1, we always negate whenever
				the column number is even (when columnNo + 1 is odd). */
				std::function<int (int) > sgn = [] (int columnNo) {
					return (columnNo % 2 == 0) ? -1 : 1;
				};

				/* Just use the minors of the entries of the top row for simplicity.
				For larger matrices, it would be more efficient to first
				find a row or column containing many zeroes to minimize the
				number of cofactor matrices. */
				for (coord_t columnNo : coords) {
					if (columnNo == 1) {
						// I know the first term uses the original sign.
						// Find bottom-right minor.
						determinant += matrix (1, columnNo)
						               * MatrixView (matrix, 2, 2, sideLength, sideLength).determinant();
					} else if (columnNo == sideLength) {
						// Find bottom-left minor.
						determinant += sgn (columnNo) * matrix (1, columnNo)
						               * MatrixView (matrix, 2, 1, sideLength, sideLength - 1).determinant();
					} else {
						// Concatenate the submatrices on either side of the current column.
						// Then take their determinant.
						determinant += sgn (columnNo) * matrix (1, columnNo) *
						               (MatrixView (matrix, 2, 1, sideLength, columnNo - 1)
						                | MatrixView (matrix, 2, columnNo + 1, sideLength, sideLength)).determinant();
					}
				}

				return determinant;
			}
		}

		//! A map of ordered pairs to values.
		numberMapType numbers;
		//! The width of the matrix, as calculated at construction/assignment.
		dimens_t width;
		//! The height of the matrix, as calculated at construction/assignment.
		dimens_t height;
	};

	Matrix::~Matrix() {}

	//Value-initialize private members!
	Matrix::Matrix() : d (new _Matrix {}) {}

	Matrix::Matrix (std::initializer_list<std::initializer_list<double>> nums)
		: Matrix()
	{
		//Iterate through every {} in nums.
		for (auto i = nums.begin(); i != nums.end(); ++i) {
			// Each {} is a row, so increment height.
			++d->height;

			/* If this is our first row, keep this row's width, and make sure
			ensuant rows have the same width. */
			if (d->width == 0) {
				d->width = i->size();
			} else if (d->width != i->size()) {
				throw std::invalid_argument ("Matrix must have uniform width.");
			}

			/*  Insert the numbers in the row. */
			for (auto j = i->begin(); j != i->end(); ++j) {
				/* Calculate the indices by subtracting iterators and adding 1
				since row and column numbers are 1-indexed. */
				d->numbers.insert (_Matrix::numberMapType::value_type
				                   ({i - nums.begin() + 1, j - i->begin() + 1}, *j));
			}
		}
	}

	Matrix::Matrix (std::initializer_list<double> nums) : Matrix ({nums})
	{
		//Construct, then transpose.
		(*this) = transposed();
	}

	Matrix::Matrix (std::istream& in) : Matrix()
	{
		// Read lines until an empty line is encountered or the stream ends.
		for (std::string line; std::getline (in, line, '\n') && !line.empty();) {
			//Each line is a new row: increment height.
			++d->height;

			std::istringstream iss;
			iss.str (line);
			decltype (d->width) numberOfColumns = 0;

			for (double a; iss >> a;) {
				++numberOfColumns;
				d->numbers.insert (_Matrix::numberMapType::value_type
				                   ({d->height, numberOfColumns}, a));
			}

			/* Only now do we know how many elements there were in the stream,
			so we compare it to a past width (if set). */
			if (d->width == 0) {
				d->width = numberOfColumns;
			} else if (d->width != numberOfColumns) {
				throw std::invalid_argument ("Matrix must have uniform width.");
			}
		}
	}

	Matrix::Matrix (const Matrix& toCopy) : Matrix()
	{
		// Just calls the copy-assignment operator
		(*this) = toCopy;
	}

	const Matrix Matrix::transposed() const
	{
		//We'll construct the transposed matrix with a stream.
		std::stringstream ss;

		//Insert by iterating by column, then by row, to exchange the coordinates.
		//i.e. the elements in one column now appear in the same row.
		for (coord_t c = 1; c <= width(); ++c) {
			for (coord_t r = 1; r <= height(); ++r) {
				ss << (*this) (r, c) << '\t';
			}

			//Insert newline to signify end of row.
			ss << std::endl;
		}

		return ss;
	}

	Matrix& Matrix::operator= (const Matrix& rhs)
	{
		//Self-assignment
		if (this == &rhs) {
			return *this;
		}

		//Copy dimensions
		d->width = rhs.width();
		d->height = rhs.height();
		//Copy numbers
		//Do not clear own map as this could invalidate MatrixView on right
		//side of assignment.
		_Matrix::numberMapType map {};

		for (coord_t r = 1; r <= rhs.height(); ++r) {
			for (coord_t c = 1; c <= rhs.width(); ++c) {
				map[ {r, c}] = rhs (r, c);
			}
		}

		d->numbers = map;
		return *this;
	}

	double& Matrix::operator() (coord_t r, coord_t c)
	{
		//Get the number at (r,c) from the map.
		return d->numbers.at ({r, c});
	}

	double Matrix::operator() (coord_t r, coord_t c) const
	{
		return d->numbers.at ({r, c});
	}

	Matrix::operator std::string() const
	{
		std::stringstream ss;

		//Insert all the numbers, separated by tabs.
		for (coord_t r = 1; r <= height(); ++r) {
			for (coord_t c = 1; c <= width(); ++c) {
				ss.width(4);
				ss.precision(3);
				if ( (*this) (r, c) == -0) {
					//Mixing negative and positive zero is ugly.
					ss << +0 << '\t';
				} else {
					ss << (*this) (r, c) << '\t';
				}
			}

			//This row is done: insert a newline.
			ss << std::endl;
		}

		return ss.str();
	}

	Matrix& Matrix::operator*= (double rhs)
	{
		for (coord_t r = 1; r <= height(); ++r) {
			for (coord_t c = 1; c <= width(); ++c) {
				//Multiply those numbers.
				(*this) (r, c) *= rhs;
			}
		}

		return *this;
	}

	const Matrix Matrix::operator* (double rhs) const
	{
		//Copy this matrix in a temporary.
		Matrix tmp = *this;
		//Perform multiplication-assignment on the temporary.
		tmp *= rhs;
		//Return the temporary.
		return tmp;
	}

	const Matrix operator* (double lhs, const Matrix& rhs)
	{
		return rhs * lhs;
	}

	Matrix& Matrix::operator/= (double rhs)
	{
		for (coord_t r = 1; r <= height(); ++r) {
			for (coord_t c = 1; c <= width(); ++c) {
				//Divide here.
				(*this) (r, c) /= rhs;
			}
		}

		return *this;
	}

	const Matrix Matrix::operator/ (double rhs) const
	{
		//Perform division-assignment on temporary copy.
		Matrix tmp = *this;
		tmp /= rhs;
		return tmp;
	};

	const Matrix Matrix::operator-() const
	{
		//Temporary copy (operator* returns temporaries) multiplied by -1.
		return (*this) * -1;
	}

	Matrix& Matrix::operator+= (const Matrix& rhs)
	{
		//You can't add them if they have different dimensions.
		if (!isSameShape (rhs)) {
			throw std::invalid_argument ("Operands must have like dimensions.");
		}

		for (coord_t r = 1; r <= height(); ++r) {
			for (coord_t c = 1; c <= width(); ++c) {
				//Add here.
				(*this) (r, c) += rhs (r, c);
			}
		}

		return *this;
	}

	const Matrix Matrix::operator+ (const Matrix& rhs) const
	{
		//Perform addition-assignment on a temporary copy.
		Matrix tmp = *this;
		tmp += rhs;
		return tmp;
	}

	Matrix& Matrix::operator-= (const Matrix& rhs)
	{
		//Negate the right hand side and then perform addition-assignment.
		*this += -rhs;
		return *this;
	}

	const Matrix Matrix::operator- (const Matrix& rhs) const
	{
		//Perform subtraction-assignment on a temporary copy.
		Matrix tmp = *this;
		tmp -= rhs;
		return tmp;
	}

	void Matrix::setAt (coord_t r, coord_t c,
	                    const Matrix& matrix)
	{
		//Iterate through given matrix
		for (coord_t rOffset = 1; rOffset <= matrix.height(); ++rOffset) {
			for (coord_t cOffset = 1; cOffset <= matrix.width(); ++cOffset) {
				//Translate to local coordinates by adding to (r,c)
				(*this) (r + rOffset - 1, c + cOffset - 1)
				    = matrix (rOffset, cOffset);
			}
		}
	}

	dimens_t Matrix::width() const
	{
		return d->width;
	}

	dimens_t Matrix::height() const
	{
		return d->height;
	}

	bool Matrix::isSameShape (const Matrix& matrix) const
	{
		//Checks if dimensions are equal.
		return width() == matrix.width() && height() == matrix.height();
	}

	const Matrix Matrix::identity (dimens_t dimension)
	{
		std::stringstream ss;

		for (dimens_t r = 1; r <= dimension; ++r) {
			//At row N, there are N - 1 zeroes before the 1
			dimens_t zeroesBeforeOne = r - 1;

			for (dimens_t i = 0; i < zeroesBeforeOne; ++i) {
				ss << 0 << '\t';
			}

			ss << 1 << '\t';

			//At row N, there are dimension - N zeroes after the 1
			dimens_t zeroesAfterOne = dimension - r;

			for (dimens_t i = 0; i < zeroesAfterOne; ++i) {
				ss << 0 << '\t';
			}

			ss << std::endl;
		}

		return ss;
	}

	Matrix& Matrix::operator*= (const Matrix& rhs)
	{
		if (width() != rhs.height()) {
			throw std::invalid_argument ("Right operand must have as many rows "
			                             "as the left has columns.");
		}

		std::stringstream ss;

		for (coord_t r = 1; r <= height(); ++r) {
			for (coord_t c = 1; c <= rhs.width(); ++c) {
				double num {};
				// A row from the left matrix
				const MatrixView leftRow = *this > "R" + std::to_string (r);
				// A column from the right matrix, transposed
				const Matrix rightTransposedColumn =
				    (rhs > "C" + std::to_string (c)).transposed();

				// Dot product
				for (coord_t i = 1; i <= width(); ++i) {
					num += leftRow (1, i) * rightTransposedColumn (1, i);
				}

				ss << num << '\t';
			}

			ss << std::endl;
		}

		Matrix tmp (ss);
		*this = tmp;
		return *this;
	}

	const Matrix Matrix::operator* (const Matrix& rhs) const
	{
		Matrix tmp = *this;
		tmp *= rhs;
		return tmp;
	}

	double Matrix::determinant() const
	{
		if (width() != height()) {
			throw std::domain_error ("Can't compute determinant of "
			                         "nonsquare matrix: "
			                         + std::to_string (width())
			                         + "x" + std::to_string (height()));
		}

		return d->determinant (*this);
	}

	const Matrix Matrix::inverse() const
	{
		if (determinant() == 0) {
			throw std::domain_error ("Matrix not invertible.");
		}
		std::cerr<<determinant() << std::endl;

		Matrix tmp = *this;
		Matrix inverse = identity (width());

		/*
		Gauss-Jordan elimination.

		For row N, divide it by the Nth element. Then, for another row M,
		subtract from row M row N multiplied by the Nth element of row M.
		An identity matrix is obtained.

		Do the same steps (but use the same scalar multiplications)
		on an identity matrix to get the inverse.
		*/
		for (coord_t r = 1; r <= height(); ++r) {
			std::string address = "R" + std::to_string (r);
			MatrixView originalRow = tmp > address;
			MatrixView inverseRow = inverse > address;
			// Divide row N by the Nth element
			double toDivide = tmp (r, r);
			originalRow /= toDivide;
			inverseRow /= toDivide;

			for (coord_t r2 = 1; r2 <= height(); ++r2) {

				std::cout << (std::string) (tmp | inverse) << std::endl <<std::endl;
				if (r2 != r) {
					std::string address2 = "R" + std::to_string (r2);
					// Subtracct from row M row N multiplied by element (M, N)
					double toMultiply = tmp (r2, r);
					(tmp > address2) -= toMultiply * originalRow;
					(inverse > address2) -= toMultiply * inverseRow;
				}
			}
		}

		return inverse;
	}

	const Matrix operator/ (const Matrix& lhs, const Matrix& rhs)
	{
		if (lhs.width() != rhs.width()) {
			throw std::invalid_argument ("Can't vertically concatenate "
			                             "matrices of different widths.");
		}

		std::stringstream ss;

		for (const Matrix & matrix : {
		lhs, rhs
	}) {
			for (dimens_t r = 1; r <= matrix.height(); ++r) {
				for (dimens_t c = 1; c <= matrix.width(); ++c) {
					// Just insert all the numbers in a row one after another.
					ss << matrix (r, c) << '\t';
				}

				ss << std::endl;
			}
		}

		return ss;
	}

	const Matrix operator| (const Matrix& lhs, const Matrix& rhs)
	{
		if (lhs.height() != rhs.height()) {
			throw std::invalid_argument ("Can't horizontally concatenate "
			                             "matrices of different heights.");
		}

		/*
		Say we have lhs:

		a b
		c d

		And rhs:

		e f
		g h

		Then lhs | rhs gives you:

		a b e f
		c d g h

		(lhs | rhs)ᵀ is:

		a c
		b d
		e g
		f h

		Which we recognize as lhsᵀ / rhsᵀ.
		*/
		return (lhs.transposed() / rhs.transposed()).transposed();
	}
};

