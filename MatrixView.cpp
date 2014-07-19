#include "Matrix.h"
#include<map>
#include<utility>
#include<memory>
#include<stdexcept>
#include<sstream>
#include<string>

namespace matrix
{
	struct MatrixView::_MatrixView {
		_MatrixView (Matrix& matrix, coord_t r1, coord_t c1,
		             coord_t r2, coord_t c2)
			: headRow (r1), headColumn (c1), width (c2 - c1 + 1),
			  height (r2 - r1 + 1), operand (matrix)
		{ }

		/*! The index of the first row in the Matrix that this view accesses. */
		coord_t headRow;
		/*! The index of the first column in the Matrix that this view accesses. */
		coord_t headColumn;
		/*! The width of the view. */
		dimens_t width;
		/*! The height of the view. */
		dimens_t height;
		/*! A pointer to a Matrix to access. */
		Matrix& operand;
	};

	MatrixView::MatrixView (const Matrix& matrix, coord_t r1, coord_t c1,
	                        coord_t r2, coord_t c2)
		: Matrix(),
		  d (new _MatrixView (const_cast<Matrix&> (matrix), r1, c1, r2, c2))
	{
		if (r1 > matrix.height() || r2 > matrix.height() ||
		        c1 > matrix.width() || c2 > matrix.width()) {
			throw std::invalid_argument ("The view's dimensions extend outside "
			                             "the matrix itself.");
		}
	}

	MatrixView::MatrixView (const MatrixView& matrixView)
		: MatrixView (matrixView.d->operand,
		              matrixView.d->headRow,
		              matrixView.d->headColumn,
		              matrixView.d->headRow + matrixView.height() - 1,
		              matrixView.d->headColumn + matrixView.width() - 1) {}

	MatrixView::~MatrixView() {}

	double& MatrixView::operator() (coord_t r, coord_t c)
	{
		// Translate to original Matrix's coordinates
		return d->operand (d->headRow + r - 1, d->headColumn + c - 1);
	}

	double MatrixView::operator() (coord_t r, coord_t c) const
	{
		return d->operand (d->headRow + r - 1, d->headColumn + c - 1);
	}

	MatrixView& MatrixView::operator= (const Matrix& rhs)
	{
		/* Since MatrixView assignment overwrites the values in another Matrix,
		it wouldn't make sense for the geometry of the new values to be different. */
		if (!isSameShape (rhs)) {
			throw std::invalid_argument ("Assignment to MatrixView requires "
			                             "matrix of same dimensions.");
		} else {
			/* Overwrite the values starting at the top-left corner of this view,
			wherever it might be on the actual Matrix. */
			setAt (1, 1, rhs);
		}

		return *this;
	}

	MatrixView& MatrixView::operator= (const MatrixView& rhs)
	{
		// Overload resolution does not work for copy assignment operators.
		return (*this) = static_cast<const Matrix&> (rhs);
	}

	dimens_t MatrixView::width() const
	{
		return d->width;
	}

	dimens_t MatrixView::height() const
	{
		return d->height;
	}

	MatrixView operator> (const Matrix& lhs, const std::string& rhs)
	{
		//See Matrix.h for the explanation of the format.
		//Extract number.
		coord_t number = std::stoll (rhs.substr (1));

		// Now do something based on the character.
		switch (rhs.at (0)) {
		case 'R':
			// A row vector
			return MatrixView (lhs, number, 1, number, lhs.width());
			break;

		case 'C':
			// A column vector
			return MatrixView (lhs, 1, number, lhs.height(), number);
			break;

		default:
			throw std::invalid_argument ("Address format incorrect: " + rhs);
		}
	}
};

