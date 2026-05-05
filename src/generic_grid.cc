#include "Interpolation/generic_grid.hh"
#include <stdexcept>
#include <algorithm>

namespace Interpolation
{
namespace Generic
{

StandardGrid::StandardGrid(const vector_d &input) : _tj(input) // passing to _tj the input
															   // vector
{
	if (input.size() <= 1) {
		throw std::invalid_argument(
			"[Generic::StandardGrid]: input vector cannot be of less than two elements.");
	}
	std::sort(_tj.begin(), _tj.end()); // sorting point

	if (_tj.front()<-1.) {
		throw std::runtime_error("Invalid vector.");
	}

	if (_tj.back()>1.) {
		throw std::runtime_error("Invalid vector.");
	}

	if (std::abs(_tj.front() - (-1.)) > 1.0e-12) {
		_tj.insert(_tj.begin(), -1.);
	} else {
		_tj.front() = -1.;
	}
	if (std::abs(_tj.back() - 1.) > 1.0e-12) {
		_tj.push_back(1.);
	} else {
		_tj.back() = 1.;
	}

	_p = _tj.size() - 1;
	_lambdaj.resize(_p+1, 1.);

	for (size_t j=0; j<=_p; j++) {
		for (size_t i=0; i<j; i++) {
			_lambdaj[j] *= (_tj[j] - _tj[i]);
		}
		for (size_t i=j+1; i<=_p; i++) {
			_lambdaj[j] *= (_tj[j] - _tj[i]);
		}
		_lambdaj[j] = 1. / _lambdaj[j];
	}

	_Dij.resize(_p+1, std::vector<double>(_p+1, 0.));

	for (size_t i=0; i<=_p; i++) {
		// diagonal element
		for (size_t n=0; n<i; n++) {
			_Dij[i][i] += 1. / (_tj[i] - _tj[n]);
		}
		for (size_t n=i+1; n<=_p; n++) {
			_Dij[i][i] += 1. / (_tj[i] - _tj[n]);
		}

		for (size_t j=0; j<=_p; j++) {
			if (j==i) continue;

			_Dij[i][j] = 1. / (_tj[i] - _tj[j]);
			for (size_t k=0; k<=_p; k++) {
				if (k==i || k==j) continue;
				_Dij[i][j] *= (_tj[j] - _tj[k]) / (_tj[i] - _tj[k]);
			}
		}
	}
} // StandardGrid::StandardGrid

StandardGrid::StandardGrid(const std::function<double(size_t, size_t)> &fnc, size_t p) : _p(p)
{
	// fnc has as a parameter the index and the total number of point
	if (std::abs(fnc(0,p)-(-1.)) > 1.0e-12) {
		throw std::invalid_argument("[GenerciGrid::StandardGrid]: lower bound must be -1");
	}
	if (std::abs(fnc(p,p)-1.) > 1.0e-12) {
		throw std::invalid_argument("[GenericGrid::StandardGrid]: upper bound must be 1");
	}

	_tj.resize(p+1, 0.);
	_tj.front() = -1.;
	_tj.back() = 1.;
	
	for (size_t i=1; i<_p; i++) {
		_tj[i] = fnc(i,p);
	}

	_lambdaj.resize(_p+1, 1.);

	for (size_t j=0; j<=_p; j++) {
		for (size_t i=0; i<j; i++) {
			_lambdaj[j] *= (_tj[j] - _tj[i]);
		}
		for (size_t i=j+1; i<=_p; i++) {
			_lambdaj[j] *= (_tj[j] - _tj[i]);
		}
		_lambdaj[j] = 1. / _lambdaj[j];
	}

	_Dij.resize(_p+1, std::vector<double>(_p+1, 0.));

	for (size_t i=0; i<=_p; i++) {
		// diagonal element
		for (size_t n=0; n<i; n++) {
			_Dij[i][i] += 1. / (_tj[i] - _tj[n]);
		}
		for (size_t n=i+1; n<=_p; n++) {
			_Dij[i][i] += 1. / (_tj[i] - _tj[n]);
		}

		for (size_t j=0; j<=_p; j++) {
			if (j==i) continue;

			_Dij[i][j] = 1. / (_tj[i] - _tj[j]);
			for (size_t k=0; k<=_p; k++) {
				if (k==i || k==j) continue;
				_Dij[i][j] *= (_tj[j] - _tj[k]) / (_tj[i] - _tj[k]);
			}
		}
	}
} // StandardGrid::StandardGrid

double StandardGrid::interpolate(double t, const vector_d &fj, size_t start, size_t end,
				STRATEGY str) const
{
	if (t<-1. || t>1. || (end-start) != _p) {
		throw std::domain_error("[GenericGrid::interpolate]: t=" + std::to_string(t)
						+ " \\notin [-1,1] OR view into fj of wrong size: ["
						+ std::to_string(start) + ", " + std::to_string(end) + "]");
	}

	switch (str) {
		case STRATEGY::NAIVE: {
			double sum = 0.;
			for (size_t i = start, j=0; i<=end; i++, j++) {
				sum += poli_weight(t, j) * fj[i];
			}			
			return sum;
			break;
		}
		case STRATEGY::FBF: {
			double monic = 1.;
			for (size_t i=0; i<=_p; i++) {
				monic *= (t - _tj[i]);
			}
			double sum = 0.;
			for (size_t i=start, j=0; i<=end; i++, j++) {
				sum += poli_weight_fbf(t, j, monic) * fj[i];
			}
			return sum;
			break;
		}
		case STRATEGY::SBF: {
			double den = 0.;
			for (size_t l=0; l<=_p; l++) {
				if (fabs(t-_tj[l]) < 1.0e-15) return fj[l+start];
				den += _lambdaj[l] / (t - _tj[l]);
			}			
			double sum = 0.;
			for (size_t i=start, j=0; i<=end; i++, j++) {
				sum += poli_weight_sbf(t, j, den) * fj[i];
			}
			return sum;
		}
        default: {
            throw std::invalid_argument("Strategy of interpolation not recognised.");
        }
	}
} // StandardGrid::interpolate
 
double StandardGrid::interpolate_der(double t, const vector_d &fj, size_t start, size_t end,
				STRATEGY str) const
{
    if (t<-1. || t>1. || (end-start) != _p) {
        throw std::domain_error("[GenericGrid::interpolate_der]: t=" + std::to_string(t)
                + " \\notin [-1,1] OR view into fj of wrong size: [" + std::to_string(start) 
                + ", " + std::to_string(end) + "]");
    }

    switch(str) {
        case STRATEGY::NAIVE: {
            double sum = 0.;
            for (size_t i=start, j=0; i<=end; i++, j++) {
                sum += poli_weight_der(t, j) * fj[i];
            }
            return sum;
            break;
        }
        case STRATEGY::FBF: {
            double monic = 1.;
            for (size_t i=0; i<=_p; i++) {
                monic *= (t - _tj[i]);
            }
            double sum = 0;
            for (size_t i=start, j=0; i<=end; i++, j++) {
                sum += poli_weight_fbf_der(t, j, monic) * fj[i];
            }
            return sum;
            break;
        }
        case STRATEGY::SBF: {
            double den = 0.;
            for (size_t l=0.; l<=_p; l++) {
                if (fabs(t-_tj[l]) < 1.0e-15) {
                    double sum = 0.;
                    for (size_t i=start, j=0; i<=end; i++, j++) {
                        sum += fj[i] * _Dij[i][l];
                    }
                    return sum;
                }
                den += _lambdaj[l] / (t - _tj[l]);
            }            
            double sum = 0.;
            for (size_t i=start, j=0; i<=end; i++, j++) {
                sum += poli_weight_sbf_der(t, j, den) * fj[i];
            }
            return sum;
        }
        default: {
            throw std::invalid_argument("Strategy of interpolation not recognised.");
        }
    }
} // StandardGrid::interpolate_der

double StandardGrid::interpolate_der_v2(double t, const vector_d &fj, size_t start, size_t end,
				STRATEGY str) const
{
    vector_d ftilde = fj;
    apply_D(ftilde, start, end);
    return interpolate(t, ftilde, start, end, str);
} // StandardGrid::interpolate_der_v2

double StandardGrid::poli_weight(double t, size_t j) const
{
    if (std::abs(t-_tj[j]) < 1.0e-15) return 1.;

    double result = 1.;
    for (size_t i=0; i<=_p; i++) {
        if (i==j) continue;
        if (std::abs(t-_tj[i]) < 1.0e-15) return 0.;
        result *= (t-_tj[i]);
    }
    return result * _lambdaj[j];
} // StandardGrid::poli_weight

double StandardGrid::poli_weight_fbf(double t, size_t j) const
{
    if (std::abs(t-_tj[j]) < 1.0e-15) return 1.;

    double monic = 1.;
    for (size_t i=0; i<=_p; i++) {
        monic *= (t - _tj[i]);
    }
    return monic * _lambdaj[j] / (t - _tj[j]);
} // StandardGrid::poli_weight_fbf

double StandardGrid::poli_weight_fbf(double t, size_t j, double monic) const
{
    if (std::abs(t-_tj[j]) < 1.0e-15) return 1.;
    return monic * _lambdaj[j] / (t - _tj[j]);
} // StandardGrid::poli_weight_fbf

double StandardGrid::poli_weight_sbf(double t, size_t j) const
{
    if (std::abs(t-_tj[j]) < 1.0e-15) return 1.;
    
    double den = 0.;
    for (size_t i=0; i<=_p; i++) {
        if (std::abs(t-_tj[i]) < 1.0e-15) return 0.;
        den += _lambdaj[i] / (t - _tj[i]);
    }
    return _lambdaj[j] / (t - _tj[j]) / den;
} // StandardGrid::poli_weight_sbf

double StandardGrid::poli_weight_sbf(double t, size_t j, double den) const
{
    if (std::abs(t-_tj[j]) < 1.0e-15) return 1.;
    return _lambdaj[j] / (t - _tj[j]) / den;
} // StandardGrid::poli_weight_sbf

double StandardGrid::poli_weight_der(double t, size_t j) const
{
    if (t<-1. || t>1.) {
        throw std::domain_error("[StandardGrid::poli_weight_der]: t=" + std::to_string(t)
                + " \\notin [-1,1]");
    }
    double res = 0.;
    for (size_t i=0; i<=_p; i++) {
        if (fabs(t - _tj[i]) < 1.0e-15) return _Dij[j][i];
        res += _Dij[j][i] * poli_weight(t, i);
    }
    return res;
} // StandardGrid::poli_weight_der

double StandardGrid::poli_weight_fbf_der(double t, size_t j) const
{
    if (t<-1 || t>1) {
        throw std::domain_error("[StandardGrid::poli_weight_sbf_der]: t=" + std::to_string(t)
                + " \\notin [-1,1]");
    }
    double res = 0.;
    for (size_t i=0; i<=_p; i++) {
        if (fabs(t - _tj[i]) < 1.0e-15) return _Dij[j][i];
        res += _Dij[j][i] * poli_weight_fbf(t, i);
    }
    return res;
} // StandardGrid::poli_weight_fbf_der
 
double StandardGrid::poli_weight_fbf_der(double t, size_t j, double monic) const
{
    if (t<-1 || t>1) {
        throw std::domain_error("[StandardGrid::poli_weight_sbf_der]: t=" + std::to_string(t)
                + " \\notin [-1,1]");
    }
    double res = 0.;
    for (size_t i=0; i<=_p; i++) {
        if (fabs(t - _tj[i]) < 1.0e-15) return _Dij[j][i];
        res += _Dij[j][i] * poli_weight_fbf(t, i, monic);
    }
    return res;
} // StandardGrid::poli_weight_fbf_der

double StandardGrid::poli_weight_sbf_der(double t, size_t j) const
{
    if (t<-1 || t>1) {
        throw std::domain_error("[StandardGrid::poli_weight_sbf_der]: t=" + std::to_string(t)
                + " \\notin [-1,1]");
    }
    double res = 0.;
    for (size_t i=0; i<=_p; i++) {
        if (fabs(t - _tj[i]) < 1.0e-15) return _Dij[j][i];
        res += _Dij[j][i] * poli_weight_sbf(t, i);
    }
    return res;
} // StandardGrid::poli_weight_sbf_der

double StandardGrid::poli_weight_sbf_der(double t, size_t j, double den) const
{
    if (t<-1 || t>1) {
        throw std::domain_error("[StandardGrid::poli_weight_sbf_der]: t=" + std::to_string(t)
                + " \\notin [-1,1]");
    }
    double res = 0.;
    for (size_t i=0; i<=_p; i++) {
        if (fabs(t - _tj[i]) < 1.0e-15) return _Dij[j][i];
        res += _Dij[j][i] * poli_weight_sbf(t, i, den);
    }
    return res;
} // StandardGrid::poli_weigth_sbf_der
 
void StandardGrid::apply_D(vector_d &fj, size_t start, size_t end) const
{
    if ((end - start) != _p) {
        throw std::invalid_argument("[StandardGrid::apply_D]: cannot apply derivative matrix"
                " to partial vector");
    }
    vector_d temp((end-start+1), 0.);
    for (size_t i=0; i<=_p; i++) {
        for (size_t j=0, k=start; k<=end; j++, k++) {
            temp[i] += fj[k] * _Dij[j][i];
        }
    }

    for (size_t i=start; i<=end; i++) {
        fj[i] = temp[i-start];
    }
} // StandardGrid::apply_D

vector_d StandardGrid::discretize(const std::function<double(double)> &fnc) const
{
    vector_d result(_p+1, 0.);
    for (size_t i=0; i<=_p; i++) {
        result[i] = fnc(_tj[i]);
    }
    return result;
} // StandardGrid::discretize

} // namespace Generic
} // namespace Interpolation
