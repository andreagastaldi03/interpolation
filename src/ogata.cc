#include "Interpolation/ogata.hh"

namespace Ogata
{

double Integrator::psi(double x)
{
   return x * std::tanh(M_PI * 0.5 * std::sinh(x));
}

double Integrator::psi_prime(double x)
{
   double num = M_PI * x * std::cosh(x) + std::sinh(M_PI * std::sinh(x));
   double den = 1. + std::cosh(M_PI * std::sinh(x));
   return num / den;
}

void Integrator::tabulate()
{
   nodes.resize(1000);
   weights.resize(1000);
   for (size_t i = 0; i < 1000; i++) {
      nodes[i] = M_PI * psi(h * J0_zero[i] / M_PI) / h;

      double tmp  = M_PI * psi_prime(h * J0_zero[i] / M_PI);
      tmp        *= std::cyl_bessel_j(0, nodes[i]);
      tmp        *= std::cyl_neumann(0, J0_zero[i]) / std::cyl_bessel_j(1, J0_zero[i]);

      weights[i] = tmp;
   }
}

double Integrator::integrate(const std::function<double(double)> &fnc, double y, double eps)
{
   double res = 0.;
   for (size_t i = 0; i < nodes.size(); i++) {
      double tmp = fnc(nodes[i] / y) * weights[i] / y;
      if (i > 10 && std::abs(tmp) / std::abs(res) < eps) break;

      res += tmp;
   }
   return res;
}

} // namespace Ogata