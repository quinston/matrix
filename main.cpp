#include<iostream>
#include<string>
#include<sstream>
#include<algorithm>
#include<cmath>
#include "Matrix.h"

std::ostream& operator<< (std::ostream& s, const matrix::Matrix& m)
{
	return s << (std::string) m;
}

int main()
{
	std::cout << "Enter one pair of X and Y values per line, the two separated "
	          "by a space: " << std::endl
	          << "(Enter a blank line when you're done.): "
	          << std::endl;

	std::stringstream ss;
	std::vector<std::pair<double, double>> xyPairs {};

	for (std::string line; std::getline (std::cin, line) && !line.empty();) {
		ss.str (line);
		double x {}, y {};
		ss >> x >> y;
		xyPairs.push_back ({x, y});
		// Clear EOF flag, since we have exhausted the stream.
		ss.clear();
	}

	ss.str ("");

	std::cout << "To what order polynomial should this data be fitted?"
	          << std::endl;
	unsigned long long order {};
	std::cin >> order;

	// Rows of X values and Y values.
	// Taking advantage of the preexisting tabularish format of Matrix's
	// string cast.
	for (std::pair<double, double> xyPair : xyPairs) {
		ss << xyPair.first << '\t' << xyPair.second << std::endl;
	}

	std::cout << std::endl << "Here is your data: " << std::endl
	          << matrix::Matrix (ss) << std::endl;
	ss.clear();
	ss.str ("");

	// The Vandermonde matrix
	for (std::pair<double, double> xyPair : xyPairs) {
		for (signed long long n = order; n >= 2; --n) {
			ss << std::pow (xyPair.first, n) << '\t';
		}

		// Don't bother with std::pow for powers 1 and 0
		ss << xyPair.first << '\t' << 1 << std::endl;
	}

	matrix::Matrix V (ss);
	std::cout << "Here is your Vandermonde matrix: " << std::endl << V
	          << std::endl;
	ss.clear();
	ss.str ("");

	std::cout << "Here is its transpose: " << std::endl << V.transposed()
	          << std::endl;

	std::cout << "Here is VᵀV: " << std::endl << V.transposed() * V
	          << std::endl;

	// Create the Y-vector
	for (std::pair<double, double> xyPair : xyPairs) {
		ss << xyPair.second << std::endl;
	}

	matrix::Matrix yVector (ss);
	ss.clear();
	ss.str ("");

	//A column vector of coefficients
	matrix::Matrix answer = (V.transposed() * V).inverse()
	                        * V.transposed()
	                        * yVector;

	std::cout << "Computed coefficients: " << std::endl;

	// Terms for power > 2
	for (matrix::coord_t r = 1; r <= answer.height() - 2; ++r) {
		std::cout << answer (r, 1) << "x^" << answer.height() - r << " + ";
	}

	// Power 1 term
	std::cout << answer (answer.height() - 1, 1) << "x + ";
	// Constant term
	std::cout << answer (answer.height(), 1) << std::endl;
	return 0;
}

