#include "Interpolation/gauss_kronrod.hh"
#include <Interpolation/grid_1d.hh>

namespace Interpolation
{

namespace
{

// u=a <=> t=+1
// u=b <=> t=-1
double from_ab_to_m1p1(double u, double a, double b) noexcept
{
   return -2 * (u - a) / (b - a) + 1;
};
double from_m1p1_to_ab(double t, double a, double b) noexcept
{
   return (b - a) * (1 - t) * 0.5 + a;
};
double from_ab_to_m1p1_der(double a, double b) noexcept
{
   return -2 / (b - a);
};

} // namespace

SingleDiscretizationInfo::SingleDiscretizationInfo(std::vector<double> inter,
                                                   std::vector<size_t> g_size,
                                                   std::function<double(double)> to_i_space,
                                                   std::function<double(double)> to_i_space_der,
                                                   std::function<double(double)> to_p_space,
                                                   std::function<double(double)> to_p_space_der)
    : intervals(inter.size() - 1, {0, 0}), intervals_phys(inter.size() - 1, {0, 0}),
      grid_sizes(g_size), to_inter_space(to_i_space), to_inter_space_der(to_i_space_der),
      to_phys_space(to_p_space), to_phys_space_der(to_p_space_der)
{
   if (g_size.size() != (inter.size() - 1)) {
      throw std::runtime_error("[SingleDiscretizationInfo] Incompatible sizes for number of "
                               "subintervals ("
                               + std::to_string(inter.size() - 1)
                               + ") andentries in the "
                                 "vector of number of points for each subinterval ("
                               + std::to_string(g_size.size()) + ")");
   }
   for (size_t i = 0; i < inter.size() - 1; i++) {
      intervals[i]      = {to_i_space(inter[i]), to_i_space(inter[i + 1])};
      intervals_phys[i] = {inter[i], inter[i + 1]};
   }
};

double Grid1D::get_der_matrix(size_t a, size_t j, size_t b, size_t k) const
{
   if (a != b) return 0.0;
   return _stored_grids.at(_d_info.grid_sizes[a])._Dij[j][k];
}

Grid1D::Grid1D(const SingleDiscretizationInfo &d_info) : _d_info(d_info)
{
   size = 0;
   if (d_info.intervals.size() == 0) {
      size_li   = 0;
      c_size    = 0;
      c_size_li = 0;
      return;
   }

   for (size_t i = 0; i < d_info.intervals.size(); i++) {
      // NOTE: This works ok because grid_sizes stores N, but the points are N+1 always
      // NOTE: This is implementation of non-overlapping grids, leading to simpler weights
      // TODO: Interpolation was good in tests, but I need to check that the Kernels are not screwed
      // by this
      size += d_info.grid_sizes[i] + 1;
      if (_stored_grids.find(d_info.grid_sizes[i]) == _stored_grids.end()) {
         _stored_grids.emplace(d_info.grid_sizes[i], Chebyshev::StandardGrid(d_info.grid_sizes[i]));
      }
   }

   size_li   = static_cast<index_t>(size);
   c_size    = size - d_info.intervals.size() + 1;
   c_size_li = static_cast<index_t>(c_size);

   _weights.resize(size);
   _weights_der.resize(size);
   _weights_sub.resize(size);
   _from_iw_to_ic.resize(size);
   _from_idx_to_inter.resize(size);

   _coord.resize(c_size);
   _coord_inter.resize(c_size);
   _delim_indexes.resize(d_info.intervals.size() + 1);

   _der_matrix.resize(size);
   for (size_t i = 0; i < size; i++)
      _der_matrix[i].resize(size);

   // No check on order of the interval
   index_t index       = 0;
   index_t index_coord = 0;
   for (size_t a = 0; a < d_info.intervals.size(); a++) {

      const Chebyshev::StandardGrid &sg = _stored_grids.at(d_info.grid_sizes[a]);
      _delim_indexes[a]                 = index;

      for (size_t j = 0; j <= d_info.grid_sizes[a]; j++) {
         // NOTE: _coord stores only the UNIQUES coordinates

         if (j != d_info.grid_sizes[a] || (a == d_info.intervals.size() - 1)) {
            _coord_inter[index_coord]
                = from_m1p1_to_ab(sg.t(j), d_info.intervals[a].first, d_info.intervals[a].second);

            _coord[index_coord] = d_info.to_phys_space(_coord_inter[index_coord]);
         }

         _from_iw_to_ic[index]     = index_coord;
         _from_idx_to_inter[index] = a;
         if (j != d_info.grid_sizes[a]) index_coord++;

         // Note: in order to have proper assignement operator
         // I cannot use [this] in the lambdas. Hence, I need to
         // pass the needed variables by value to the lambdas.
         size_t cached_int_size                 = _d_info.intervals.size();
         std::pair<double, double> cached_inter = _d_info.intervals[a];

         // ------------------------------
         _weights[index] = [a, j, cached_int_size,
                            cached_inter](double u, const Chebyshev::StandardGrid &sg) -> double {
            double res = 0;

            bool condition_x
                = a == cached_int_size - 1 ? u <= cached_inter.second : u < cached_inter.second;
            if (u >= cached_inter.first && condition_x) {
               res += sg.poli_weight(from_ab_to_m1p1(u, cached_inter.first, cached_inter.second),
                                     j);
            }

            return res;
         };

         // ------------------------------
         _weights_der[index] = [a, j, cached_int_size, cached_inter](
                                   double u, const Chebyshev::StandardGrid &sg) -> double {
            double res = 0;

            bool condition_x
                = a == cached_int_size - 1 ? u <= cached_inter.second : u < cached_inter.second;
            if (u >= cached_inter.first && condition_x) {
               const double dl_dx = from_ab_to_m1p1_der(cached_inter.first, cached_inter.second);
               const double dw_dl = sg.poli_weight_der(
                   from_ab_to_m1p1(u, cached_inter.first, cached_inter.second), j);
               res += dw_dl * dl_dx;
            }

            return res;
         };
         _weights_sub[index] = [a, j, cached_int_size, cached_inter](
                                   double u, const Chebyshev::StandardGrid &sg) -> double {
            double res = 0;

            bool condition_x
                = a == cached_int_size - 1 ? u <= cached_inter.second : u < cached_inter.second;
            if (u >= cached_inter.first && condition_x) {
               res += sg.poli_weight(from_ab_to_m1p1(u, cached_inter.first, cached_inter.second), j)
                    - 1;
            }

            return res;
         };

         index_t inner_index = 0;
         for (size_t b = 0; b < d_info.intervals.size(); b++) {
            for (size_t k = 0; k <= d_info.grid_sizes[b]; k++) {
               _der_matrix[index][inner_index++] = get_der_matrix(a, j, b, k);
            }
         }
         index++;
      }
   }
   _delim_indexes[d_info.intervals.size()] = index;

   _from_ic_to_iw.resize(c_size);
   for (size_t i = 0; i < _from_iw_to_ic.size(); i++) {
      _from_ic_to_iw[_from_iw_to_ic[i]].push_back(i);
   }

   {
      using integrator = GaussKronrod<GK_41>;

      _integral_weights.resize(c_size_li, 0.);

      size_t i_w;

      std::function<double(double)> full_integrand = [&](double v) -> double {
         return _d_info.to_phys_space_der(v) * _weights[i_w](v, get_std_grid(i_w));
      };

      for (size_t j = 0; j < size; j++) {

         size_t j_c = _from_iw_to_ic[j];

         i_w = j;

         auto [vmin, vmax] = get_support_weight_aj(j);

         if (std::fabs(vmax - vmin) < 1.0e-15) {
            continue;
         }
         _integral_weights[j_c]
             += integrator::integrate(full_integrand, vmin, vmax, 1.0e-10, 1.0e-10);

      } // j-loop
   }
}

} // namespace Interpolation