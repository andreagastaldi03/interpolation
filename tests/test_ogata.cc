#include <Interpolation/interpolation.hh>
#include <iostream>
int main()
{
   {
      auto fnc = [](double b) {
         return b * exp(-b * b / 2.) / 2.;
      };

      Ogata::Integrator integrator;
      std::cout << integrator.integrate(fnc, 1., 1.0e-10) - 0.3032653298563167 << std::endl;
      std::cout << integrator.integrate(fnc, 2., 1.0e-10) - 0.06766764161830635 << std::endl;
   }
}