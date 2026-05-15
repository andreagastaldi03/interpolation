#include <Interpolation/interpolation.hh>
#include <iostream>

using namespace Interpolation;

double foo(double x)
{
   return log(x);
}

int main()
{
   Grid1D grid(make_discretization_info<details::log_0_maps>({1.0e-3, 0.2, 1.0}, {16, 16}));
   std::vector<double> fj = Discretize<std::vector<double>, double>(grid, foo, [](size_t n) {
      return std::vector<double>(n, 0.);
   });

   double y = 0.56661151;

   double tmp = grid.interpolate<double, std::vector<double>>(y, fj, []() -> double {
      return 0.;
   });

   double exact = foo(y);

   std::printf("%.16e\n", tmp - exact);
   std::printf("%.16e\n", tmp);
   std::printf("%.16e\n", exact);

   tmp = 0.;
   std::cout << "--" << std::endl;

   for (size_t i = 0; i < fj.size(); i++) {
      tmp += fj[i] * grid._integral_weights[i];
      std::cout << fj[i] << " " << grid._integral_weights[i] << std::endl;
   }
   std::cout << "--" << std::endl;

   exact = -0.9920922447210179;

   std::printf("%.16e\n", tmp - exact);
   std::printf("%.16e\n", tmp);
   std::printf("%.16e\n", exact);
}