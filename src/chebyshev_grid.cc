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
		if (std::abs(t - _tj[j]) < 1.e-15) return fj[j + start]; //if t is a grid point
		// you already know the value of the approximation (=f evaluated in that point),
		// so thank to that return you close directly the function interpolate, not 
		// only the for.
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
{
	if (t < -1 || t > 1 || (end-start) != _p) {
		throw std::domain_error("[StandardGrid::interpolate]: t=" + std::to_string(t)
					+ " \\notin [-1,1] OR view into fj of wrong size: ["
					+ std::to_string(start) + ", " + std::to_string(end) + "]");
	}

	double den = 0.;
	for (size_t l=0; l<=_p; l++) {
		if (std::abs(t-_tj[l])<1.0e-15) {
			double sum = 0;
			for (size_t i = start, j=0; i <= end; i++, j++) {
				sum += fj[i] * _Dij[j][l];
			}
			return sum;
		}
		den += _betaj[l] / (t - _tj[l]);
	}

	double sum = 0.;
	for (size_t i = start, j=0; i<=end; i++, j++) {
		sum += poli_weight_der(t, j, den) * fj[i];
		// it's writing the function f' as sum over i of fj l'j, and l'j is calculated
		// in poli_weight_der
	}
	return sum;
} // StandardGrid::interpolate_der
 
double StandardGrid::poli_weight(double t, size_t j) const
{
	if (t < -1 || t > 1) {
      throw std::domain_error("[StandardGrid::poli_weight]: t=" + std::to_string(t)
                              + " \\notin [-1, +1]");
	}

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
{
	if (t < -1 || t > 1) {
		throw std::domain_error("[StandardGrid::poli_weight]: t=" + std::to_string(t)
                              + " \\notin [-1, +1]");
	}

	double res = 0;
	// it's calculating l'j(x)
	for (size_t i = 0; i<=_p; i++) {
		if (std::abs(t - _tj[i]) < 1.0e-15) return _Dij[j][i]; // l'j(xi) is Dji
		res += _Dij[j][i] * poli_weight(t,i); // if not a grid point, write l'j as
		// a series using lj, l'j sum over i of l'j(xi) li = Dji li
	}
	return res;
} // StandardGrid::poli_weight_der
 
double StandardGrid::poli_weight_der(double t, size_t j, double den) const
{
	if (t < -1 || t > 1) {
		throw std::domain_error("[StandardGrid::poli_weight]: t=" + std::to_string(t)
                              + " \\notin [-1, +1]");
	}

	double res = 0;
	for (size_t i = 0; i<=_p; i++) {
		if (std::abs(t - _tj[i]) < 1.0e-15) return _Dij[j][i]; 
		res += _Dij[j][i] * poli_weight(t,i,den); 
	}
	return res;
} // StandardGrid::poli_weight_der

void StandardGrid::apply_D(vector_d &fj, size_t start, size_t end) const
{
	if (end - start != _p) {
      throw std::invalid_argument("[StandardGrid::apply_D]: cannot apply "
                                  "derivative matrix to partial vector.");
    }

	vector_d temp((end-start+1), 0.);
	for (size_t i=0; i<= _p; i++) {
		for (size_t j=0, k=start; k<=end; j++, k++) {
			temp[i] += fj[k] * _Dij[j][i]; // f'(xi) = sum over j fj l'j(xi)
		}
	}

	for (size_t i=start; i<=end; i++) {
		fj[i] = temp[i-start];
	}
} // StandardGrid::apply_D

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
