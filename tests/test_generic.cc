#include "Interpolation/interpolation.hh"
using namespace Interpolation;

double testfunction(double x)
{
   return exp(2 * x);
}

double testfunction_d(double x)
{
   return 2 * exp(2 * x);
}

double equi_pts(size_t j, size_t p)
{
   return 2 * ((double)j / ((double)p)) - 1;
}

int test_exp_function()
{
   using ST = Generic::StandardGrid::STRATEGY;

   const auto testfunction = [](double x) {
      return exp(2.0 * x);
   };
   const auto testfunction_d = [](double x) {
      return 2.0 * exp(2.0 * x);
   };

   const double p = 20;
   Generic::StandardGrid grid(equi_pts, p);

   int errcode = 1;

   if (grid._p != p) return errcode;
   if (grid._lambdaj.size() != p + 1) return errcode;
   if (grid._tj.size() != p + 1) return errcode;
   if (grid._Dij.size() != p + 1) return errcode;
   for (size_t i = 0; i <= p; i++) {
      if (grid._Dij[i].size() != p + 1) return errcode;
   }

   auto v = grid.discretize(testfunction);

   double xmin = -1;
   double xmax = 1;
   size_t n    = 1e4;
   double dx   = (xmax - xmin) / ((double)n - 1);

   double m_n = 0, md_n = 0.;
   double m_1 = 0, md_1 = 0.;
   double m_2 = 0, md_2 = 0.;

   std::FILE *fptr = std::fopen("StandardGrid_interpolation_chebyshev.dat", "w");
   for (size_t i = 0; i < n; i++) {
      const double x       = xmin + i * dx;
      const double exact   = testfunction(x);
      const double inter_n = grid.interpolate(x, v, 0, grid._p, ST::NAIVE);
      const double inter_1 = grid.interpolate(x, v, 0, grid._p, ST::FBF);
      const double inter_2 = grid.interpolate(x, v, 0, grid._p, ST::SBF);

      const double exact_d   = testfunction_d(x);
      const double inter_d_n = grid.interpolate_der_v2(x, v, 0, grid._p, ST::NAIVE);
      const double inter_d_1 = grid.interpolate_der_v2(x, v, 0, grid._p, ST::FBF);
      const double inter_d_2 = grid.interpolate_der_v2(x, v, 0, grid._p, ST::SBF);

      m_n  = std::max(m_n, std::abs(exact - inter_n));
      md_n = std::max(md_n, std::abs(exact_d - inter_d_n));

      m_1  = std::max(m_1, std::abs(exact - inter_1));
      md_1 = std::max(md_1, std::abs(exact_d - inter_d_1));

      m_2  = std::max(m_2, std::abs(exact - inter_2));
      md_2 = std::max(md_2, std::abs(exact_d - inter_d_2));

      std::fprintf(fptr, "%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\n", x,
                   exact, inter_n, inter_1, inter_2, exact_d, inter_d_n, inter_d_1, inter_d_2);
   }
   for (size_t i = 0; i <= grid._p; i++) {
      const double x       = grid._tj[i];
      const double exact   = testfunction(x);
      const double inter_n = grid.interpolate(x, v, 0, grid._p, ST::NAIVE);
      const double inter_1 = grid.interpolate(x, v, 0, grid._p, ST::FBF);
      const double inter_2 = grid.interpolate(x, v, 0, grid._p, ST::SBF);

      const double exact_d   = testfunction_d(x);
      const double inter_d_n = grid.interpolate_der_v2(x, v, 0, grid._p, ST::NAIVE);
      const double inter_d_1 = grid.interpolate_der_v2(x, v, 0, grid._p, ST::FBF);
      const double inter_d_2 = grid.interpolate_der_v2(x, v, 0, grid._p, ST::SBF);

      m_n  = std::max(m_n, std::abs(exact - inter_n));
      md_n = std::max(md_n, std::abs(exact_d - inter_d_n));

      m_1  = std::max(m_1, std::abs(exact - inter_1));
      md_1 = std::max(md_1, std::abs(exact_d - inter_d_1));

      m_2  = std::max(m_2, std::abs(exact - inter_2));
      md_2 = std::max(md_2, std::abs(exact_d - inter_d_2));
   }
   std::fclose(fptr);

   std::printf("Naive Max difference in interpolation: %.6e\n", m_n);
   std::printf("Naive Max difference in interpolation derivatives: %.6e\n", md_n);

   std::printf("FBF Max difference in interpolation: %.6e\n", m_1);
   std::printf("FBF Max difference in interpolation derivatives: %.6e\n", md_1);

   std::printf("SBF Max difference in interpolation: %.6e\n", m_2);
   std::printf("SBF Max difference in interpolation derivatives: %.6e\n", md_2);

   if (m_n > 1.0e-11) return errcode;
   if (md_n > 1.0e-10) return errcode;

   if (m_1 > 1.0e-11) return errcode;
   if (md_1 > 1.0e-10) return errcode;

   if (m_2 > 1.0e-11) return errcode;
   if (md_2 > 1.0e-10) return errcode;

   return 0;
}

int main()
{
   int i;
   i = test_exp_function();
   if (i != 0) return i;
   return 0;
}
