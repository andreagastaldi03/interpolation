#include "Interpolation/interpolation.hh"
#include <iostream>

struct DYGrid {
   DYGrid(size_t nOgata, size_t nXi) : nO(nOgata), nXi(nXi), data(nOgata * nXi, 0.)
   {
   }

   double &operator()(size_t i, size_t j)
   {
      return data[i * nXi + j];
   }

   const double &operator()(size_t i, size_t j) const
   {
      return data[i * nXi + j];
   }

   double *block(size_t i)
   {
      size_t offset = i * nXi;
      return data.data() + offset;
   }

   size_t nO, nXi;
   std::vector<double> data;
};

struct Table {
   std::vector<double> qT;
   std::vector<double> ogata_grid;
   std::vector<double> xi_grid;
   std::vector<DYGrid> WDY;
};

Table compute(const std::vector<double> &qT, std::pair<double, double> bin, size_t nOgata)
{
   using namespace Interpolation;

   if (bin.first >= bin.second) {
      throw std::logic_error("bin.first cannot be >= bin.second");
   }

   std::vector<double> bounds;
   std::vector<size_t> sizes;
   if (bin.first < 1 && bin.second > 1) {
      bounds = {bin.first, 1., bin.second};
      sizes  = {10, 10};
   } else {
      bounds = {bin.first, bin.second};
      sizes  = {20};
   }

   Grid1D grid(make_discretization_info<details::log_0_maps>(bounds, sizes));

   using integrator = GaussKronrod<GK_21>;
   Ogata::Integrator ogata_integrator;

   Table table;
   table.xi_grid = grid._coord;
   size_t nXi    = table.xi_grid.size();
   table.qT      = qT;
   for (size_t i = 0; i < nOgata; i++) {
      table.ogata_grid.push_back(ogata_integrator.nodes[i]);
   }

   for (size_t i = 0; i < qT.size(); i++) {
      table.WDY.emplace_back(DYGrid(nOgata, nXi));

      for (size_t iO = 0; iO < nOgata; iO++) {
         double *dygrid = table.WDY.back().block(iO);

         const double b = ogata_integrator.nodes[iO] / qT[i];
         int i_xi       = 0;

         const double pref = b * ogata_integrator.weights[iO] / (0.1 + qT[i] * qT[i]);

         auto integrand = [&i_xi, &grid, b](double u) {
            const double xi  = grid._d_info.to_phys_space(u);
            const double jac = grid._d_info.to_phys_space_der(u);
            const double w   = grid._weights[i_xi](u, grid.get_std_grid(i_xi));

            return exp(-b * b) * pow(xi, -0.25 /*-1+0.75*/) * w * jac;
         };

         for (size_t jXi = 0; jXi < grid.size; jXi++) {
            i_xi = jXi;

            auto [umin, umax] = grid.get_support_weight_aj(jXi);

            dygrid[grid._from_iw_to_ic[jXi]] += pref * integrator::integrate(integrand, umin, umax);
         }
      }
   }

   return table;
}

std::vector<double> convolve(const Table &tab, const std::function<double(double, double)> &fnc)
{
   std::vector<double> result(tab.qT.size(), 0.);

   for (size_t iqT = 0; iqT < tab.qT.size(); iqT++) {
      const double qT = tab.qT[iqT];

      double &r            = result[iqT];
      const DYGrid &dygrid = tab.WDY[iqT];

      for (size_t iO = 0; iO < tab.ogata_grid.size(); iO++) {
         const double b = tab.ogata_grid[iO] / qT;

         double acc = 0.;

         for (size_t iXi = 0; iXi < tab.xi_grid.size(); iXi++) {
            acc += dygrid(iO, iXi) * fnc(tab.xi_grid[iXi], b) * fnc(1. / tab.xi_grid[iXi], b);
         }

         r += acc;
         if ((std::abs(acc / r) < 1.0e-5 && iO > 10) || (std::abs(r) < 1e-16 && iO > 10)) break;
      }
   }

   return result;
}

int main()
{

   {
      auto fnc = [](double xi, double b) {
         return xi * (1 - xi) / (1. + b * b);
      };

      const std::vector<double> qT = {0.1, 0.2, 0.3, 0.4, 0.5};

      auto table = compute(qT, {0.1, 2.2}, 200);

      auto result = convolve(table, fnc);

      std::vector<double> exacts = {-0.03520884293764395, -0.07016612787804205,
                                    -0.10462372668327038, -0.138340298216793, -0.17108450219757113};

      for (size_t i = 0; i < result.size(); i++) {
         std::printf("%.2f\t%.6e\n", table.qT[i], result[i] - exacts[i]);
      }
   }

   {
      auto fnc = [](double xi, double b) {
         return xi * (0.01 - xi) / (1. + b * b);
      };

      std::vector<double> qT;
      for (size_t i = 1; i < 20; i++) {
         qT.push_back((double)i / 20.);
      }

      auto table = compute(qT, {1., exp(2.4)}, 200);

      auto result = convolve(table, fnc);

      std::FILE *fp = std::fopen("CS.dat", "w");
      for (size_t i = 0; i < result.size(); i++) {
         std::fprintf(fp, "%.16e\t%.16e\n", table.qT[i], result[i]);
      }
      std::fclose(fp);
   }

   if (false) {
      using namespace Interpolation;
      Grid1D grid(make_discretization_info<details::log_0_maps>({0.1, 1, 2.2}, {10, 10}));

      size_t i_xi;
      auto integrand = [&i_xi, &grid](double u) {
         const double xi  = grid._d_info.to_phys_space(u);
         const double jac = grid._d_info.to_phys_space_der(u);
         const double w   = grid._weights[i_xi](u, grid.get_std_grid(i_xi));

         return pow(xi, -0.25 /*-1+0.75*/) * w * jac;
      };
      using integrator = GaussKronrod<GK_21>;

      double res = 0.;

      auto fnc = [](double xi) {
         return xi * (1 - xi);
      };

      for (size_t jXi = 0; jXi < grid.size; jXi++) {
         i_xi = jXi;

         auto [umin, umax] = grid.get_support_weight_aj(jXi);

         const double xi  = grid._coord[grid._from_iw_to_ic[jXi]];
         res             += integrator::integrate(integrand, umin, umax) * fnc(xi) * fnc(1. / xi);
      }
      std::cout << res << std::endl;
   }

   return 0;
}