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

constexpr double _LOCAL_DBL_MIN_     = 2.2250738585072014e-308;
constexpr double _1e8_LOCAL_DBL_MIN_ = 2.2250738585072014e-300;
constexpr double _LOCAL_DBL_EPSILON_ = 2.2204460492503131e-16;
constexpr double _LOCAL_DBL_MAX_     = 1.7976931348623157e+308;

static inline bool subinterval_too_small(double a1, double a2, double b2)
{
   const double tmp = (1 + _LOCAL_DBL_EPSILON_) * (std::fabs(a2) + _LOCAL_DBL_MIN_);
   return std::fabs(a1) <= tmp && std::fabs(b2) <= tmp;
}

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

} // namespace Interpolation