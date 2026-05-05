#pragma once

#include "Interpolation/default.hh"
#include <cmath>

namespace Interpolation
{
namespace details
{

/**
 * @brief Implements the map \f$ u = \log\left(\log\left(\frac{\mu^2}{\Lambda^2}\rigth)\right) \f$
 * with \f$ \Lambda^2 = 0.2^2 \f$. The implementation assumes that the input physical variable is
 * \f$ \log(\mu^2) = t \f$ directly
 *
 */
struct log_log_mu2_maps {
   /// Log(Lambda^2) = Log(0.2^2)
   static constexpr double LLambda2 = -3.2188758248682006;

   /// To inteprolation space
   static double tis(double t)
   {
      return std::log(t - LLambda2);
   }
   /// To interpolation space derivative
   static double tis_d(double t)
   {
      return 1.0 / (t - LLambda2);
   }
   /// To physical space
   static double tps(double z)
   {
      return std::exp(z) + LLambda2;
   }
   /// To physical space derivative
   static double tps_d(double z)
   {
      return std::exp(z);
   }
};

/// Identity map \f$ u = z \f$
struct identity_maps {

   /// To inteprolation space
   static double tis(double t)
   {
      return t;
   }
   /// To inteprolation space derivative
   static double tis_d(double)
   {
      return 1.0;
   }
   /// To physical space
   static double tps(double z)
   {
      return z;
   }
   /// To physical space derivative
   static double tps_d(double)
   {
      return 1.0;
   }
};

/// Logarithmic map towards lower limit, with regularized end points
struct log_0_maps {
   static constexpr double eps = 1.0e-3;

   /// To inteprolation space
   static double tis(double t)
   {
      return log((t + eps) / (1 + eps));
   }
   /// To inteprolation space derivative
   static double tis_d(double t)
   {
      return 1.0 / (t + eps);
   }
   /// To physical space
   static double tps(double z)
   {
      return exp(z) * (1 + eps) - eps;
   }
   /// To physical space derivative
   static double tps_d(double z)
   {
      return exp(z) * (1 + eps);
   }
};

/// Logarithmic map towards upper limit, with regularized end points
struct log_1_maps {
   static constexpr double eps = 1.0e-3;

   /// To inteprolation space
   static double tis(double t)
   {
      return -std::log((1.0 - eps) * (1 - t) + eps);
   }
   /// To inteprolation space derivative
   static double tis_d(double t)
   {
      return (1.0 - eps) / ((1.0 - eps) * (1 - t) + eps);
   }
   /// To physical space
   static double tps(double z)
   {
      return (1.0 - std::exp(-z)) / (1.0 - eps);
   }
   /// To physical space derivative
   static double tps_d(double z)
   {
      return std::exp(-z) / (1.0 - eps);
   }
};

struct atanh_maps {
   static constexpr double eps = 1.0e-3;

   /// To inteprolation space
   static double tis(double t)
   {
      return atanh(-1. + 2. * ((1.0 - 2. * eps) * t + eps));
   }
   /// To inteprolation space derivative
   static double tis_d(double t)
   {
      return 2.0 * (1.0 - 2. * eps) / (1.0 - pow((1. - 2. * t) * (1.0 - 2. * eps), 2));
   }
   /// To physical space
   static double tps(double z)
   {
      return (1. - 2. * eps + tanh(z)) / (2. - 4 * eps);
   }
   /// To physical space derivative
   static double tps_d(double z)
   {
      return 1.0 / pow(cosh(z), 2) / (2. - 4. * eps);
   }
};

template <double P, double E>
struct powlaw_0_maps {
   static constexpr double p   = P;
   static constexpr double eps = E;

   /// To inteprolation space
   static double tis(double t)
   {
      return -pow(t + eps, p);
   }
   /// To inteprolation space derivative
   static double tis_d(double t)
   {
      return -p * pow(t + eps, p - 1.);
   }
   /// To physical space
   static double tps(double z)
   {
      return pow(-z, 1. / p) - eps;
   }
   /// To physical space derivative
   static double tps_d(double z)
   {
      return -pow(-z, 1. / p - 1.) / p;
   }
};

} // namespace details

} // namespace Interpolation
