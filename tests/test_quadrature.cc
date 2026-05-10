#include <Interpolation/gauss_kronrod.hh>
#include <cmath>
#include <iostream>
#include <iomanip>

void integate_and_print(const std::function<double(double)> &fnc, const std::string &name,
                        const double exact)
{
   using integrator = Interpolation::GaussKronrod<Interpolation::GK_21>;

   const double res1 = integrator::integrate(fnc, 0, 1, 1.0e-10, 1.0e-10);
   const double res2 = integrator::integrate_rec(fnc, 0, 1, 1.0e-10, 1.0e-10);

   std::cout << std::setprecision(6) << std::scientific;
   std::cout << name << "\t";
   std::cout << "Global adaptive: " << res1 - exact << "\t";
   std::cout << "Local adaptive: " << res2 - exact << "\t";
   std::cout << std::endl;
}

int main()
{
   {
      auto fnc = [](double x) {
         return exp(x);
      };
      integate_and_print(fnc, "exp", exp(1) - 1);
   }
   {
      for (size_t i = 0; i < 50; i += 5) {
         auto fnc = [i](double x) {
            return std::pow(x, i);
         };
         integate_and_print(fnc, "x^" + std::to_string(i), 1. / (i + 1.));
      }
   }

   {
      auto fnc = [](double x) {
         return std::cyl_bessel_j(0, 10 * x);
      };
      integate_and_print(fnc, "J_0(10 x) ", 0.10670113039567368);
   }
}