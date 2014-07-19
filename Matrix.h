#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include<iosfwd>
#include<string>
#include<iterator>
#include<initializer_list>
#include<memory>

/*! Includes the matrix class and related functions and classes. */
namespace matrix
{
	//! A coordinate type.
	typedef unsigned long long coord_t;
	//! A dimension (e.g. width) type.
	typedef unsigned long long dimens_t;
	/*! Matrix class. Does matrix math. Cannot be resized except by
	reassignment.
	*/
	class Matrix
	{
	public:
		/*! Constructs an empty matrix. */
		Matrix();
		/*! Constructs from a list of rows. */
		Matrix (std::initializer_list<std::initializer_list<double>>);
		/*! Constructs a column vector from a list of numbers. */
		Matrix (std::initializer_list<double>);
		/*! Continuously reads from the istream until a newline is found, which
		determines the number of columns. Keeps reading until the end of stream.
		Throws an exception if there aren't enough numbers in subsequent rows. */
		Matrix (std::istream&);
		/*! Copies the values in the given matrix. */
		Matrix (const Matrix&);
		/*! Destructor. */
		virtual ~Matrix();

		/*! Calculates the transpose.
		\return A copy of the matrix, transposed. */
		const Matrix transposed() const;
		/*! Calculates the inverse.
		Descriptive exceptions are thrown if the matrix isn't invertible.
		\return A copy of the matrix inverted. */
		const Matrix inverse() const;
		/*! Generates an identity matrix with the given dimension. */
		static const Matrix identity (dimens_t);
		/*! Calculates the determinant. */
		double determinant() const;

		/*! Matrix multiplication. */
		const Matrix operator* (const Matrix&) const;
		/*! Multiplication-assignment, overloaded for matrix multiplication. */
		Matrix& operator*= (const Matrix&);

		/*! Scalar multiplication. */
		const Matrix operator* (double) const;
		/*! Multiplication-assignment, for a scalar. */
		Matrix& operator*= (double);
		/*! Scalar division. */
		const Matrix operator/ (double) const;
		/*! Division-assignment, for a scalar. */
		Matrix& operator/= (double);
		/*! Negation of all components. */
		const Matrix operator-() const;

		/*! Adds all the components.
		Throws an exception if the matrices have different shapes. */
		const Matrix operator+ (const Matrix&) const;
		/*! Addition-assignment. */
		Matrix& operator+= (const Matrix&);
		/*! Subtracts all the components.
		Throws an exception if the matrices have different shapes. */
		const Matrix operator- (const Matrix&) const;
		/*! Subtraction-assignment. */
		Matrix& operator-= (const Matrix&);

		/*! Copy assignment. Replaces previous values completely. */
		virtual Matrix& operator= (const Matrix&);

		/*! Starting at row \em r, column \em c, fill in the given matrix's
		values. If the calling matrix is too small, throws an exception. */
		void setAt (coord_t r, coord_t c, const Matrix&);

		/*! Returns a string representing the matrix. */
		operator std::string() const;

		/*! Gets the element at the \em rth row's \em cth column. */
		virtual double& operator() (coord_t r, coord_t c);

		/*! \copydoc Matrix::operator()(coord_t,coord_t) */
		virtual double operator() (coord_t r, coord_t c) const;

		/*! Gets the number of rows. */
		virtual dimens_t height() const;

		/*! Gets the number of columns. */
		virtual dimens_t width() const;

		/*! Checks if this and another Matrix have the same dimensions. */
		bool isSameShape (const Matrix&) const;

	private:
		//! Private implementation class
		class _Matrix;
		std::unique_ptr<_Matrix> d;
	};

	/*! A view on some subset of a Matrix. When modified, the modifications
	are reflected in the original Matrix, and vice versa if the original
	Matrix is modified. The target area of a MatrixView cannot be changed. */
	class MatrixView : public Matrix
	{
	public:
		/*! Constructs a view from row \em r1, column \em c1 of the given
		matrix to the row \em r2, column \em c2. If the dimensions are
		out of bounds, an exception is thrown. */
		MatrixView (const Matrix&, coord_t r1, coord_t c1, coord_t r2,
		            coord_t c2);
		/*! Copy constructor. */
		MatrixView (const MatrixView&);
		/*! Destructor. */
		virtual ~MatrixView();

		/*! Gets the element at this view's \em rth row's \em cth column. */
		double& operator() (coord_t r, coord_t c);

		/*! \copydoc MatrixView::operator()(coord_t,coord_t) */
		double operator() (coord_t r, coord_t c) const;

		/*! Replaces the view's values with the given matrix. If the sizes are
		different, an exception is thrown. */
		MatrixView& operator= (const Matrix&);
		/*! \copydoc MatrixView::operator=(const Matrix&)

		This is explicitly defined since T = T (T being the type) always
		attempt to call operator=(T) and will not accept operator=(A) where
		A is some base of T. */
		MatrixView& operator= (const MatrixView&);

		/*! \copydoc Matrix::width() */
		dimens_t width() const;

		/*! \copydoc Matrix::height() */
		dimens_t height() const;

	private:
		//! Private implementation class
		class _MatrixView;
		std::unique_ptr<_MatrixView> d;
	};

	//! \copydoc Matrix::operator*(double) const
	const Matrix operator* (double, const Matrix&);

	/*! Returns a MatrixView given an address of the following format:

	\em thing \em number

	Where \em thing is one of the following:
	\par R
	Row
	\par C
	Column

	And \em number is the row or column number (one-indexed).
	*/
	MatrixView operator> (const Matrix&, const std::string&);

	/*! Vertically concatenates 2 matrices of the same width. */
	const Matrix operator/ (const Matrix&, const Matrix&);

	/*! Horizontally concatenates 2 matrices of the same height.  */
	const Matrix operator| (const Matrix&, const Matrix&);
};
#endif // MATRIX_H_INCLUDED

