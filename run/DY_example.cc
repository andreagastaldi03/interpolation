#include "Interpolation/interpolation.hh"
#include <iostream>
#include "Interpolation/ran2.hh"
#include <ceres/ceres.h>

Ran2 rng(-2);

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

std::pair<std::vector<double>, std::vector<std::vector<double>>>
convolve_der(const Table &tab, const std::function<double(double, double)> &fnc,
             const std::function<void(double, double, std::vector<double> &)> &fnc_der,
             size_t n_par)
{
   std::vector<std::vector<double>> result_der(tab.qT.size(), std::vector<double>(n_par, 0.));
   std::vector<double> result(tab.qT.size(), 0.);

   std::vector<double> cache_1(n_par, 0.);
   std::vector<double> cache_2(n_par, 0.);

   for (size_t iqT = 0; iqT < tab.qT.size(); iqT++) {
      const double qT = tab.qT[iqT];

      auto &r_der          = result_der[iqT];
      double &r            = result[iqT];
      const DYGrid &dygrid = tab.WDY[iqT];
      std::vector<int> should_stop(n_par + 1, false);

      for (size_t iO = 0; iO < tab.ogata_grid.size(); iO++) {
         const double b = tab.ogata_grid[iO] / qT;

         double acc = 0.;
         std::vector<double> acc_der(n_par, 0.);

         for (size_t iXi = 0; iXi < tab.xi_grid.size(); iXi++) {
            fnc_der(tab.xi_grid[iXi], b, cache_1);
            fnc_der(1. / tab.xi_grid[iXi], b, cache_2);

            double f1 = fnc(tab.xi_grid[iXi], b);
            double f2 = fnc(1. / tab.xi_grid[iXi], b);

            acc += dygrid(iO, iXi) * f1 * f2;
            for (size_t iPar = 0; iPar < n_par; iPar++) {
               acc_der[iPar] += dygrid(iO, iXi) * (cache_1[iPar] * f2 + f1 * cache_2[iPar]);
            }
         }

         int stop_count = 0;
         for (size_t iPar = 0; iPar < n_par; iPar++) {
            if (!should_stop[iPar]) {
               r_der[iPar]       += acc_der[iPar];
               should_stop[iPar]  = ((std::abs(acc_der[iPar] / r_der[iPar]) < 1.0e-5 && iO > 10)
                                    || (std::abs(r_der[iPar]) < 1e-16 && iO > 10));
            } else stop_count++;
         }
         if (!should_stop[n_par]) {
            r += acc;
            should_stop[n_par]
                = (std::abs(acc / r) < 1.0e-5 && iO > 10) || (std::abs(r) < 1e-16 && iO > 10);
         } else stop_count++;

         if ((size_t)stop_count == n_par + 1) break;
      }
   }

   return {result, result_der};
}

const std::vector<double> qT_generated = {
    5.0000000000000003e-02, 1.0000000000000001e-01, 1.4999999999999999e-01, 2.0000000000000001e-01,
    2.5000000000000000e-01, 2.9999999999999999e-01, 3.4999999999999998e-01, 4.0000000000000002e-01,
    4.5000000000000001e-01, 5.0000000000000000e-01, 5.5000000000000004e-01, 5.9999999999999998e-01,
    6.5000000000000002e-01, 6.9999999999999996e-01, 7.5000000000000000e-01, 8.0000000000000004e-01,
    8.4999999999999998e-01, 9.0000000000000002e-01, 9.4999999999999996e-01};

const std::vector<double> data_generated = {
    6.2387823210035165e-01, 1.1618903278774217e+00, 1.5626645503759611e+00, 1.8193070765924098e+00,
    1.9540043957417932e+00, 1.9988632325181352e+00, 1.9836820620421454e+00, 1.9314426878536923e+00,
    1.8581745338909843e+00, 1.7743901691132440e+00, 1.6866914516515581e+00, 1.5990771812254057e+00,
    1.5138605525936442e+00, 1.4322995282435893e+00, 1.3550028175236111e+00, 1.2821776030597598e+00,
    1.2138021472230098e+00, 1.1497190458304070e+00, 1.0897028174637247e+00};

std::vector<double> fluctuate_data()
{
   std::vector<double> data(data_generated.size(), 0.0);
   for (size_t i = 0; i < data_generated.size(); i++) {
      data[i] = data_generated[i] + rng.normal(0, 0.1);
   }
   return data;
}

struct Chi2 {
   using model_t = std::function<double(double, double, const std::vector<double> &)>;
   using model_der_t
       = std::function<void(double, double, const std::vector<double> &, std::vector<double> &)>;

   Chi2(const model_t &model, const model_der_t &model_der,
        const std::vector<double> &initial_params)
       : model(model), model_der(model_der), data(fluctuate_data()), params(initial_params)
   {
      table = compute(qT_generated, {1., exp(2.4)}, 200);
   }

   void FillResidual(double *residual)
   {
      const auto pred = convolve(table, [this](double xi, double b) {
         return model(xi, b, params);
      });

      for (size_t i = 0; i < pred.size(); i++) {
         residual[i] = (pred[i] - data[i]) / std::sqrt(0.01);
      }
   }

   void FillResidualAndDerive(double *residual, double **jacobians)
   {
      const auto [pred, pred_der] = convolve_der(
          table,
          [this](double xi, double b) {
             return model(xi, b, params);
          },
          [this](double xi, double b, std::vector<double> &out) {
             return model_der(xi, b, params, out);
          },
          params.size());
      for (size_t i = 0; i < pred.size(); i++) {
         residual[i] = (pred[i] - data[i]) / std::sqrt(0.01);
         for (size_t j = 0; j < params.size(); j++) {
            jacobians[j][i] = pred_der[i][j] / std::sqrt(0.01);
         }
      }
   }

   double Evaluate()
   {
      const auto pred = convolve(table, [this](double xi, double b) {
         return model(xi, b, params);
      });
      double res      = 0.;
      for (size_t i = 0; i < pred.size(); i++) {
         res += (pred[i] - data[i]) * (pred[i] - data[i]) / 0.01;
      }
      return res;
   }

   std::vector<double> Derive()
   {
      const auto [pred, pred_der] = convolve_der(
          table,
          [this](double xi, double b) {
             return model(xi, b, params);
          },
          [this](double xi, double b, std::vector<double> &out) {
             return model_der(xi, b, params, out);
          },
          params.size());

      std::vector<double> res(params.size(), 0.);
      for (size_t i = 0; i < pred.size(); i++) {
         for (size_t j = 0; j < params.size(); j++) {
            res[j] += 2. * pred_der[i][j] * (pred[i] - data[i]) / 0.01;
         }
      }
      return res;
   }

   model_t model;
   model_der_t model_der;
   std::vector<double> data;
   std::vector<double> params;
   Table table;
};

class CeresLossFunction : public ceres::CostFunction
{
public:
   CeresLossFunction(Chi2 &chi2) : _chi2(chi2)
   {
      set_num_residuals(chi2.data.size());
      const size_t npar = chi2.params.size();
      _pbs.assign(npar, 1);

      _ceres_params.reserve(_pbs.size());
      for (const size_t &block_size : _pbs) {
         mutable_parameter_block_sizes()->push_back(static_cast<int32_t>(block_size));
         _ceres_params.push_back(new double[block_size]);
      }

      size_t curr_ip = 0;
      curr_ip        = 0;
      for (size_t i = 0; i < _pbs.size(); i++) {
         for (size_t j = 0; j < _pbs[i]; j++) {
            _index_map[curr_ip] = std::make_pair(i, j);
            _ceres_params[i][j] = chi2.params[curr_ip++];
         }
      }
   }

   bool Evaluate(double const *const *parameters, double *residuals,
                 double **jacobians) const override
   {
      for (size_t i = 0; i < _chi2.get().params.size(); i++) {
         std::pair<size_t, size_t> rc = _index_map.at(i);
         _chi2.get().params[i]        = parameters[rc.first][rc.second];
      }

      if (jacobians != nullptr) {
         _chi2.get().FillResidualAndDerive(residuals, jacobians);
      } else {
         _chi2.get().FillResidual(residuals);
      }
      return true;
   }

   /// To make C++ happy
   ~CeresLossFunction() override {};

   std::vector<double> GetCurrentParameters() const
   {
      std::vector<double> params(_chi2.get().params.size(), 0.);
      for (size_t i = 0; i < params.size(); i++) {
         std::pair<size_t, size_t> rc = _index_map.at(i);

         params[i] = _ceres_params[rc.first][rc.second];
      }
      return params;
   }

public:
   std::reference_wrapper<Chi2> _chi2;
   /// Vector of parameters, ceres required format
   std::vector<double *> _ceres_params;
   /// Parameter block sizes
   std::vector<size_t> _pbs;
   /// Index map, from flatten parameter index to pair of [block_index, position_in_block]
   std::map<size_t, std::pair<size_t, size_t>> _index_map;
};

int main()
{

   std::vector<std::vector<double>> min_params;

   std::vector<double> loss_values;

   for (size_t n_repl = 0; n_repl < 1000; n_repl++) {
      auto fnc = [](double xi, double b, const std::vector<double> &params) {
         return xi * (params[0] - xi) / (1. + params[1] * b * b);
      };
      auto fnc_der
          = [](double xi, double b, const std::vector<double> &params, std::vector<double> &out) {
               double v = xi * (params[0] - xi) / (1. + params[1] * b * b);
               out[0]   = xi / (1. + params[1] * b * b);
               out[1]   = -v * b * b / (1. + params[1] * b * b);
            };

      Chi2 chi2(fnc, fnc_der, {0.5, 0.5});

      CeresLossFunction loss(chi2);
      std::unique_ptr<ceres::Problem> _problem;
      {
         ceres::Problem::Options problem_options;
         problem_options.cost_function_ownership = ceres::DO_NOT_TAKE_OWNERSHIP;
         _problem = std::make_unique<ceres::Problem>(problem_options);
         _problem->AddResidualBlock(&loss, nullptr, loss._ceres_params);
      }
      ceres::Solver::Options _options;
      {
         _options.max_num_iterations           = 100;
         _options.minimizer_progress_to_stdout = false;
         _options.function_tolerance           = 1.0e-15;
         _options.gradient_tolerance           = 1.0e-6;
         _options.parameter_tolerance          = 1.0e-6;
      }
      ceres::Solver::Summary _summary;

      ceres::Solve(_options, _problem.get(), &_summary);

      auto final_params = loss.GetCurrentParameters();
      min_params.push_back(final_params);
      chi2.params = final_params;
      loss_values.emplace_back(chi2.Evaluate());

      std::cout << n_repl << std::endl;
   }

   std::FILE *fp;

   {
      fp = std::fopen("Parameters.dat", "w");
      for (size_t i = 0; i < min_params.size(); i++) {
         for (size_t j = 0; j < min_params[i].size(); j++) {
            std::fprintf(fp, "%.16e\t", min_params[i][j]);
         }
         std::fprintf(fp, "\n");
      }
      std::fclose(fp);
   }

   {
      fp = std::fopen("Losses.dat", "w");
      for (size_t i = 0; i < loss_values.size(); i++) {
         std::fprintf(fp, "%.16e\n", loss_values[i]);
      }
      std::fclose(fp);
   }
}

int main_tests()
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