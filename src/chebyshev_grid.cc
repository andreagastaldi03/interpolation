#include "Interpolation/chebyshev_grid.hh"
#include <stdexcept>

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

double StandardGrid::interpolate(double t, const vector_d &fj, size_t start, size_t end) const
{
	if (t < -1 || t > 1) {
		throw std::domain_error("StandardGrid::interpolate t must be in [-1,1]"); 
		// calling smt outside the domain 
	}
	if (end - start != _p) {
		throw std::domain_error("StandardGrid::interpolate end-start should be = p");
	}
	double den = 0.;
	for (size_t j=0; j<=_p; j++) {
		// if (t == _tj[j]) return fj[j+start];
		if (std::abs(t - _tj[j]) < 1.e-15) return fj[j + start];
		den += _betaj[j] / (t - _tj[j]);
	}	
	
	double res = 0.;
	for (size_t i=0; i<=_p; i++) {
		res += poli_weight(t, i, den) * fj[i+start];
		// res += poli_weight(t, i) * fj[i+start];
	}
	return res;
} // StandardGrid::interpolate

double StandardGrid::interpolate_der(double t, const vector_d &fj, size_t start, size_t end) const
{} // StandardGrid::interpolate_der
 
double StandardGrid::poli_weight(double t, size_t j) const
{
	if (std::abs(t - _tj[j]) < 1.e-15) return 1.;

	double den = 0.;
	for (size_t j=0; j<=_p; j++) {
		if (std::abs(t - _tj[j]) < 1.e-15) return 0.;
		den += _betaj[j] / (t - _tj[j]);
	}	

	double res = 0.;
	res        = _betaj[j] / (t - _tj[j]) / den;
	
	return res;
} // StandardGrid::poli_weight
 
double StandardGrid::poli_weight(double t, size_t j, double den) const
{
	if (std::abs(t - _tj[j]) < 1.e-15) return 1.;

	double res = 0.;
	res        = _betaj[j] / (t - _tj[j]) / den;
	
	return res;
} // StandardGrid::poli_weight
 
double StandardGrid::poli_weight_der(double t, size_t j) const
{} // StandardGrid::poli_weight_der
 
double StandardGrid::poli_weight_der(double t, size_t j, double den) const
{} // StandardGrid::poli_weight_der

void StandardGrid::apply_D(vector_d &fj, size_t start, size_t end) const
{} // StandardGrid::apply_D

vector_d StandardGrid::discretize(const std::function<double(double)> &fnc) const
{
	vector_d fj(_p+1, 0.);
	for (size_t i=0; i<=_p; i++) {
		fj[i] = fnc(_tj[i]);
	}

	return fj;
} // StandardGrid::discretize

} // namespace Chebyshev
} // namespace Interpolation
