#include "Interpolation/interpolation.hh"
#include <iostream>

double testfnc(double x)
{
	return exp(x);
}

int main()
{
	Interpolation::Chebyshev::StandardGrid grid(25);

	auto fj = grid.discretize(testfnc);

	double res = grid.interpolate(0.27, fj, 0, fj.size() - 1);
	double exact = testfnc(0.27);

	std::cout << res << std::endl;
	std::cout << exact << std::endl;
	std::cout << res - exact << std::endl;

	return 0;
}
