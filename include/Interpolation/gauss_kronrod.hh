#pragma once

#include "Interpolation/default.hh"
#include <queue>

namespace Interpolation
{

namespace cpt
{

template <typename Rule, std::size_t N>
concept IsGKRule = requires(Rule r) {
   { r.size } -> std::convertible_to<size_t>;
   requires std::same_as<decltype(r.x_gk), const std::array<double, N>>;
   requires std::same_as<decltype(r.w_g), const std::array<double, (N - 1) / 2>>;
   requires std::same_as<decltype(r.w_gk), const std::array<double, N>>;
};

} // namespace cpt

struct GK_21 {
   constexpr static size_t size = 11;
   const static std::array<double, 11> x_gk;
   const static std::array<double, 5> w_g;
   const static std::array<double, 11> w_gk;
};

struct GK_41 {
   constexpr static size_t size = 21;
   const static std::array<double, 21> x_gk;
   const static std::array<double, 10> w_g;
   const static std::array<double, 21> w_gk;
};

struct GK_61 {
   constexpr static size_t size = 31;
   const static std::array<double, 31> x_gk;
   const static std::array<double, 15> w_g;
   const static std::array<double, 31> w_gk;
};

template <typename Rule>
requires cpt::IsGKRule<Rule, Rule::size>
struct GaussKronrod {

   using rule_t = Rule;

   /**
    * @brief Integrate the input function over an interval
    *
    * @param fnc The function to be integrated
    * @param a The lower limit of the integral
    * @param b The upper limit of the integral
    * @param toll_rel The maximum relative error accepted
    * @param toll_abs The maximum absolute error accepted
    * @return double
    */
   static double integrate(const std::function<double(double)> &fnc, double a, double b,
                           double toll_rel = 1.0e-6, double toll_abs = 1.0e-10);

   static double integrate_rec(const std::function<double(double)> &fnc, double a, double b,
                               double toll_rel = 1.0e-6, double toll_abs = 1.0e-10);

private:
   struct Eval {
      double res, err;
   };
   struct Item {
      Eval e;
      double low, high;
   };
   struct ByErrAbsMax {
      bool operator()(const Item &lhs, const Item &rhs) const
      {
         return lhs.e.err < rhs.e.err;
      }
   };
   using PQ = std::priority_queue<Item, std::vector<Item>, ByErrAbsMax>;

   static Eval gauss_kronrod_simplified(const std::function<double(double)> &fnc, double a,
                                        double b);
};

namespace details
{

constexpr double _LOCAL_DBL_MIN_     = 2.2250738585072014e-308;
constexpr double _1e8_LOCAL_DBL_MIN_ = 2.2250738585072014e-300;
constexpr double _LOCAL_DBL_EPSILON_ = 2.2204460492503131e-16;
constexpr double _LOCAL_DBL_MAX_     = 1.7976931348623157e+308;

static inline bool subinterval_too_small(double a1, double a2, double b2)
{
   const double tmp = (1 + _LOCAL_DBL_EPSILON_) * (std::fabs(a2) + _LOCAL_DBL_MIN_);
   return std::fabs(a1) <= tmp && std::fabs(b2) <= tmp;
}

static inline double rescale_error(double err, const double result_abs, const double result_asc)
{
   if (result_asc != 0. && err != 0.) {
      double t = (200 * err / result_asc);
      err      = result_asc * std::min(t * std::sqrt(t), 1.);
   }
   if (result_abs > _LOCAL_DBL_MIN_ / (50 * _LOCAL_DBL_EPSILON_)) {
      double min_err = 50 * _LOCAL_DBL_EPSILON_ * result_abs;
      if (min_err > err) return min_err;
      else return err;
   }

   return err;
}
} // namespace details

template <typename Rule>
requires cpt::IsGKRule<Rule, Rule::size>
auto GaussKronrod<Rule>::gauss_kronrod_simplified(const std::function<double(double)> &fnc,
                                                  double a, double b) -> GaussKronrod<Rule>::Eval
{
   const size_t n               = Rule::size;
   const double center          = 0.5 * (a + b);
   const double half_length     = 0.5 * (b - a);
   const double abs_half_length = std::fabs(half_length);
   const double f_center        = fnc(center);

   std::array<double, Rule::size> a_fv1{};
   std::array<double, Rule::size> a_fv2{};

   double result_gauss   = 0;
   double result_kronrod = f_center * Rule::w_gk[n - 1];
   double result_abs     = std::fabs(result_kronrod);
   double result_asc     = 0.;

   if constexpr (Rule::size % 2 == 0) result_gauss = f_center * Rule::w_g[n / 2 - 1];

   for (size_t j = 0, jtw = 1; j < (n - 1) / 2; j++, jtw += 2) {
      const double abscissa  = half_length * Rule::x_gk[jtw];
      const double fv1       = fnc(center - abscissa);
      const double fv2       = fnc(center + abscissa);
      const double fsum      = fv1 + fv2;
      a_fv1[jtw]             = fv1;
      a_fv2[jtw]             = fv2;
      result_gauss          += Rule::w_g[j] * fsum;
      result_kronrod        += Rule::w_gk[jtw] * fsum;
      result_abs            += Rule::w_gk[jtw] * (std::fabs(fv1) + std::fabs(fv2));
   }

   for (size_t j = 0; j < n - 1; j += 2) {
      const double abscissa  = half_length * Rule::x_gk[j];
      const double fv1       = fnc(center - abscissa);
      const double fv2       = fnc(center + abscissa);
      a_fv1[j]               = fv1;
      a_fv2[j]               = fv2;
      result_kronrod        += Rule::w_gk[j] * (fv1 + fv2);
      result_abs            += Rule::w_gk[j] * (std::fabs(fv1) + std::fabs(fv2));
   }

   double mean = result_kronrod * 0.5;
   result_asc  = Rule::w_gk[n - 1] * std::fabs(f_center - mean);

   for (size_t j = 0; j < n - 1; j++) {
      result_asc += Rule::w_gk[j] * (std::fabs(a_fv1[j] - mean) + std::fabs(a_fv2[j] - mean));
   }

   double err = std::fabs(result_kronrod - result_gauss) * half_length;

   result_kronrod *= half_length;
   result_abs     *= abs_half_length;
   result_asc     *= abs_half_length;

   return Eval{.res = result_kronrod, .err = details::rescale_error(err, result_abs, result_asc)};
}

template <typename Rule>
requires cpt::IsGKRule<Rule, Rule::size>
double GaussKronrod<Rule>::integrate(const std::function<double(double)> &fnc, double A, double B,
                                     double toll_rel, double toll_abs)
{

   /// Hard coded for the moment
   std::size_t max_refine    = 200000;
   std::size_t max_intervals = 300000;
   PQ pq;

   double curr_res = 0.0;
   double curr_err = 0.0;

   size_t too_small_cnt = 0;

   // seed with the whole interval
   {
      Eval e = gauss_kronrod_simplified(fnc, A, B);
      pq.emplace(Item{.e = e, .low = A, .high = B});
      curr_res += e.res;
      curr_err += e.err;
   }

   for (size_t iters = 0; iters < max_refine; iters++) {
      if (pq.empty()) break;
      if (pq.size() > max_intervals) {
         // Add warning message
         break;
      }
      if (!std::isfinite(curr_res)) {
         // Add error message
         break;
      }
      if (!std::isfinite(curr_err)) {
         // Add error message
         break;
      }

      Item current = pq.top();
      pq.pop();

      double a1, b1, a2, b2;

      a1 = current.low;
      b1 = 0.5 * (current.low + current.high);
      a2 = b1;
      b2 = current.high;

      if (details::subinterval_too_small(a1, b1, b2)) {
         // Do not push back the interval, so that is no longer in the queue for
         // further, useless, refinements. Do not update the total result so that
         // this interval has its contribution 'frozen'
         if ((++too_small_cnt) >= 10) {
            // Add warning message
         }
         continue;
      }

      Eval first  = gauss_kronrod_simplified(fnc, a1, b1);
      Eval second = gauss_kronrod_simplified(fnc, a2, b2);

      pq.emplace(Item{.e = first, .low = a1, .high = b1});
      pq.emplace(Item{.e = second, .low = a2, .high = b2});

      curr_res += first.res + second.res - current.e.res;
      curr_err += first.err + second.err - current.e.err;

      double tolerance = std::max(toll_abs, toll_rel * std::fabs(curr_res));

      if (curr_err <= tolerance) break;
   }

   return curr_res;
}

template <typename Rule>
requires cpt::IsGKRule<Rule, Rule::size>
double GaussKronrod<Rule>::integrate_rec(const std::function<double(double)> &fnc, double A,
                                         double B, double toll_rel, double toll_abs)
{

   struct Inter {
      Eval e;
      double low, high;
      size_t depth;
   };

   /// Hard coded for the moment

   /// This is already quite deep recursion
   std::size_t max_depth = 50;

   std::size_t max_refine    = 200000;
   std::size_t max_intervals = 300000;

   std::vector<Inter> stack;

   double curr_res = 0.0;

   // seed with the whole interval
   {
      Eval e = gauss_kronrod_simplified(fnc, A, B);
      stack.emplace_back(Inter{.e = e, .low = A, .high = B, .depth = 0});
      curr_res += e.res;
   }

   for (size_t iters = 0; iters < max_refine; iters++) {
      if (stack.empty()) break;
      if (stack.size() > max_intervals) {
         // Add warning message
         break;
      }
      if (!std::isfinite(curr_res)) {
         // Add error message
         break;
      }

      Inter current = stack.back();
      stack.pop_back();

      double a1, b1, a2, b2;

      a1 = current.low;
      b1 = 0.5 * (current.low + current.high);
      a2 = b1;
      b2 = current.high;

      if (details::subinterval_too_small(a1, b1, b2)) {
         // Do not push back the interval, so that is no longer in the queue for
         // further, useless, refinements. Do not update the total result so that
         // this interval has its contribution 'frozen'
         continue;
      }

      if (current.depth > max_depth) {
         // We subdivided too many times a single interval.
         continue;
      }

      Eval first  = gauss_kronrod_simplified(fnc, a1, b1);
      Eval second = gauss_kronrod_simplified(fnc, a2, b2);

      double tolerance_1 = std::max(toll_abs, toll_rel * std::fabs(first.res));
      double tolerance_2 = std::max(toll_abs, toll_rel * std::fabs(second.res));

      // Push back in the stack only if tollerance is not reached.
      if (first.err > tolerance_1) {
         stack.emplace_back(Inter{.e = first, .low = a1, .high = b1, .depth = current.depth});
      }
      if (second.err > tolerance_2) {
         stack.emplace_back(Inter{.e = second, .low = a2, .high = b2, .depth = current.depth});
      }

      curr_res += first.res + second.res - current.e.res;
   }

   return curr_res;
}

} // namespace Interpolation

