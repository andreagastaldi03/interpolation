#include "Interpolation/chebyshev_grid.hh"

namespace Interpolation
{
namespace Chebyshev
{
StandardGrid::StandardGrid(size_t p)
{
	 // this is a constructor, 
	 // so i have to define everything 
	 // in the struct
	
	_p = p;

	// t_j = cos(j pi / p) for j = 0, .. , p
	for (size_t j = 0; j <= p; j++) {
		_tj.push_back(cos(j * M_PI / static_cast<double>(p)));
	}

	// beta_j = (-1)^j * (1 if j!=0, p otherwise 0.5)
	for (size_t j = 0; j <= p; j++) {
	    double sign = j % 2 == 0 ? +1 : -1; // if construction
									   // is j divisible by 2?
									   // yes -> +1, otherwise -1
		double scaling = 1.;
		if (j==0 || j==p) scaling = 0.5; // rewriting 
		_betaj.push_back(sign * scaling);
	}

	_Dij.resize(p+1, vector_d(p+1, 0.));
	// D is a vector of dim p+1 OF vector of dim p+1 OF zeros
	_Dij[0][0] = (2. * p * p + 1.) / 6.;
	_Dij[p][p] = - _Dij[0][0];

	for (size_t j = 1; j < p; j++) {
		_Dij[j][j] = - 0.5 * _tj[j] / (1. - _tj[j] * _tj[j]);
	}

	for (size_t i = 0; i <= p; i++) {
		for (size_t j = 0; j <= p; j++) {
			if (j == i) continue;
			_Dij[i][j] = -(_betaj[i] / _betaj[j]) / (_tj[i] - _tj[j]);
		}	
	}
} // StandardGrid::StandardGrid(size_t p)

} // namespace Chebyshev
} // namespace Interpolation
