#include "Interpolation/interpolation.hh"
#include <iostream>

int main()
{
	Interpolation::Chebyshev::StandardGrid grid(4); // constructor

	for (size_t j=0; j < grid._betaj.size(); j++) {
		std::cout << grid._betaj[j] << std::endl;
	}

	return 0;
}
