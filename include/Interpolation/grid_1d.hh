#pragma once

#include "Interpolation/default.hh"
#include "Interpolation/chebyshev_grid.hh"
#include "Interpolation/generic_grid.hh"
#include "Interpolation/grid_maps.hh"
#include "Interpolation/grid_concepts.hh"
#include <cstdint>

namespace Interpolation
{
struct SingleDiscretizationInfo {
public:
   /// @name Default Constructors
   /// @{
   /// @cond
   SingleDiscretizationInfo()                                                     = default;
   SingleDiscretizationInfo(const SingleDiscretizationInfo &other)                = default;
   SingleDiscretizationInfo &operator=(const SingleDiscretizationInfo &other)     = default;
   SingleDiscretizationInfo(SingleDiscretizationInfo &&other) noexcept            = default;
   SingleDiscretizationInfo &operator=(SingleDiscretizationInfo &&other) noexcept = default;
   /// @endcond
   ///@}

   /**
    * @brief Construct a new SingleDiscretizationInfo from the specified parameters
    *
    * @param inter           Intervals division points in the physical variable
    * @param g_size          Polynomial degree in each interval
    * @param to_i_space      Map from physical space to interpolation space (default: identity)
    * @param to_i_space_der  Derivative map from physical space to interpolation space (default:
    * 'one' fnc)
    * @param to_p_space      Map from interpolation space to physical space (default: identity)
    * @param to_p_space_der  Derivative map from interpolation space to physical space (default:
    * 'one' fnc)
    */
   SingleDiscretizationInfo(std::vector<double> inter, std::vector<size_t> g_size,
                            std::function<double(double)> to_i_space = details::identity_maps::tis,
                            std::function<double(double)> to_i_space_der
                            = details::identity_maps::tis_d,
                            std::function<double(double)> to_p_space = details::identity_maps::tps,
                            std::function<double(double)> to_p_space_der
                            = details::identity_maps::tps_d);

   /// Intervals in interpolation space
   std::vector<std::pair<double, double>> intervals;
   /// Intervals in physical space
   std::vector<std::pair<double, double>> intervals_phys;
   /// Polynomial degree for each interval
   std::vector<size_t> grid_sizes;

   /// Map to interpolation space
   std::function<double(double)> to_inter_space;
   /// Derivative of map to interpolation space
   std::function<double(double)> to_inter_space_der;
   /// Map to physical space
   std::function<double(double)> to_phys_space;
   /// Derivative of map to physical space
   std::function<double(double)> to_phys_space_der;
};

template <class PIM>
requires cpt::isPIM<PIM>
static SingleDiscretizationInfo make_discretization_info(std::vector<double> inter,
                                                         std::vector<size_t> g_size)
{
   return SingleDiscretizationInfo(inter, g_size, PIM::tis, PIM::tis_d, PIM::tps, PIM::tps_d);
};

/**
 * @brief Grid1D struct
 *
 * Containes all the grid points and weights that characterize a grid
 */
struct Grid1D {
   /// @name Default Constructors
   /// @{
   /// Compiler supplied constructors
   /// @cond
   Grid1D()                                   = default;
   Grid1D(const Grid1D &other)                = default;
   Grid1D &operator=(const Grid1D &other)     = default;
   Grid1D(Grid1D &&other) noexcept            = default;
   Grid1D &operator=(Grid1D &&other) noexcept = default;
   /// @endcond
   /// @}

   /**
    * @brief Construct a new Grid1D object
    *
    * @param d_info The Discretization info structure, holds the intervals, the polynomial
    * degrees and the physical <-> interpolation space maps
    */
   Grid1D(const SingleDiscretizationInfo &d_info);

   /**
    * @brief Get the der matrix entry D_{a,j; b,k}
    *
    * @param a left sub-interval index
    * @param j index in the a-th sub-interval
    * @param b right sub-interval index
    * @param k index in the b-th sub-interval
    * @return double D_{a,j; b,k} (= 0 if a!=b)
    */
   double get_der_matrix(size_t a, size_t j, size_t b, size_t k) const;

   /**
    * @brief Get the support of the `index`-th weight in interpolation space
    *
    * @param index The index of the weight
    * @return std::pair<double, double> The support of the weight (low, high)
    */
   std::pair<double, double> get_support_weight_aj(size_t index) const
   {
      return _d_info.intervals[_from_idx_to_inter[index]];
   }

   /**
    * @brief Get the support of the `index`-th weight in physical space
    *
    * @param index The index of the weight
    * @return std::pair<double, double> The support of the weight (low, high)
    */
   std::pair<double, double> get_phys_support_weight_aj(size_t index) const
   {
      return _d_info.intervals_phys[_from_idx_to_inter[index]];
   }

   /**
    * @brief Get Reference to standard grid for the given weight index
    *
    */
   const Chebyshev::StandardGrid &get_std_grid(size_t j) const
   {
      return _stored_grids.at(_d_info.grid_sizes[_from_idx_to_inter[j]]);
   }

   /**
    * @brief Interpolate on the grid
    *
    * @tparam ReturnType The type of the output variable
    * @tparam InnerType The type of the inner variable (indexable container, with elements of
    * type ReturnType)
    * @param y the point in which the interpolation is requested
    * @param input The input discretization container, contains the discretized values on the
    * grid
    * @param make_zero Create a `zero` object of type ReturnType
    * @return The interpolated value
    *
    * @warning ReturnType and InnerType must satisfy the `cpt::InterpolateCompatible<ReturnType,
    * InnerType>` concept, which ensures that the InnerType can be properly indexed, multiplied
    * by scalar values, and its element can be summed to a ReturnType variable
    */
   template <class ReturnType, class InnerType, class... Args>
   requires cpt::InterpolateCompatible<ReturnType, InnerType, Args...>
   ReturnType interpolate(double y, const InnerType &input,
                          const std::function<ReturnType()> &make_zero, Args &&...args) const
   {
      const double u = _d_info.to_inter_space(y);
      ReturnType res = make_zero();
      for (size_t j = 0; j < size; j++) {
         auto supp = get_support_weight_aj(j);
         if (u < supp.first || u > supp.second) continue;
         if constexpr (cpt::InterpolateCompatibleIndex<ReturnType, InnerType>) {
            res += input[_from_iw_to_ic[j]] * _weights[j](u, get_std_grid(j));
         } else if constexpr (cpt::InterpolateCompatibleEvaluate<ReturnType, InnerType, Args...>) {
            res += input[_from_iw_to_ic[j]].Evaluate(std::forward<Args>(args)...)
                 * _weights[j](u, get_std_grid(j));
         } else if constexpr (cpt::InterpolateCompatibleFunctor<ReturnType, InnerType, Args...>) {
            res += input[_from_iw_to_ic[j]](std::forward<Args>(args)...)
                 * _weights[j](u, get_std_grid(j));
         } else {
            throw std::invalid_argument("Cannot interpolate");
         }
      }
      return res;
   }

   /**
    * @brief Interpolate the derivative on the grid
    *
    * @tparam ReturnType The type of the output variable
    * @tparam InnerType The type of the inner variable (indexable container, with elements of
    * type ReturnType)
    * @param y The point in which the interpolation is requested
    * @param jac_ext The external jacobian, for any additional variable transformation
    * @param input The input discretization container, contains the discretized values on the
    * grid
    * @param make_zero Create a `zero` object of type ReturnType
    * @return The derivative of the interpolated value
    *
    * @warning ReturnType and InnerType must satisfy the `cpt::InterpolateCompatible<ReturnType,
    * InnerType>` concept, which ensures that the InnerType can be properly indexed, multiplied
    * by scalar values, and its element can be summed to a ReturnType variable
    */
   template <class ReturnType, class InnerType, class... Args>
   requires cpt::InterpolateCompatible<ReturnType, InnerType, Args...>
   ReturnType interpolate_der(double y, const InnerType &input,
                              const std::function<ReturnType()> &make_zero, Args &&...args) const
   {
      const double u   = _d_info.to_inter_space(y);
      const double jac = _d_info.to_inter_space_der(y);
      ReturnType res   = make_zero();
      for (size_t j = 0; j < size; j++) {
         auto supp = get_support_weight_aj(j);
         if (u < supp.first || u > supp.second) continue;
         if constexpr (cpt::InterpolateCompatibleIndex<ReturnType, InnerType>) {
            res += input[_from_iw_to_ic[j]] * jac * _weights_der[j](u, get_std_grid(j));
         } else if constexpr (cpt::InterpolateCompatibleEvaluate<ReturnType, InnerType, Args...>) {
            res += input[_from_iw_to_ic[j]].Evaluate(std::forward<Args>(args)...) * jac
                 * _weights_der[j](u, get_std_grid(j));
         } else if constexpr (cpt::InterpolateCompatibleFunctor<ReturnType, InnerType, Args...>) {
            res += input[_from_iw_to_ic[j]](std::forward<Args>(args)...) * jac
                 * _weights_der[j](u, get_std_grid(j));
         } else {
            throw std::invalid_argument("Cannot interpolate");
         }
      }
      return res;
   }

   /**
    * @brief Get the double weights object
    *
    * @return Eigen::MatrixXd int w_i w_j
    */
   matrix_d get_double_weights();

   /**
    * @brief Computes the discretized integral on the grid
    *
    * @param integrand The integrand function, depends on (x, y)
    * @param weight_fnc Eventual weight Q function of (y)
    * @return Eigen::MatrixXd The discretized representation of
    * the integral
    * \f[
    *   \int_l^u dy K(x_i, y) Q(y) w_j(y) \equiv M_{ij}
    * \f]
    */
   matrix_d integrate(
       const std::function<double(double, double)> &integrand,
       const std::function<double(double)> &weight_fnc = [](double) {
          return 1.;
       }) const;

   /**
    * @brief Computes the integral on the grid
    *
    * @param integrand The integrand function, depends on (y)
    * @return vector_d The discretized representation of
    * the integral
    * \f[
    *   \int_l^u dy K(y) w_j(y) \equiv V_{j}
    * \f]
    */
   vector_d integrate(const std::function<double(double)> &integrand, double eps = 1.0e-10) const;

   /**
    * @brief Computes the integral on the grid, subtracted at the lower bound
    *
    * @param integrand The integrand function, depends on (y)
    * @return vector_d The discretized representation of
    * the integral
    * \f[
    *   \int_l^u dy K(y) w_j(y) \equiv V_{j}
    * \f]
    * for \f$ j\neq 0 \f$ and
    * \f[
    *   \left[\int_l^u dy K(y) (w_0(y)-1)\right] + \int_l^u dy K(y)
    * \f]
    */
   vector_d integrate_subtr(const std::function<double(double)> &integrand) const;

   /// Discretization infos
   SingleDiscretizationInfo _d_info;
   /// Stored StandardGrids
   std::map<size_t, Chebyshev::StandardGrid> _stored_grids;
   /// The intepolating weights
   std::vector<std::function<double(double, const Chebyshev::StandardGrid &)>> _weights;
   /// Derivative of the interpolating weights (dw / du)
   std::vector<std::function<double(double, const Chebyshev::StandardGrid &)>> _weights_der;
   /// Weights minus 1 (w-1)
   std::vector<std::function<double(double, const Chebyshev::StandardGrid &)>> _weights_sub;
   /// Grid1D nodes in physical space
   std::vector<double> _coord;
   /// Grid1D nodes in interpolation space
   std::vector<double> _coord_inter;
   /// Conversion map from index in the weight array to index in the node array
   std::vector<i32> _from_iw_to_ic;
   /// Conversion map from index in the node array to indexes in the weight array
   std::vector<std::vector<i32>> _from_ic_to_iw;
   /// Conversion map from index in the weight array to sub-interval index (idx -> a)
   std::vector<size_t> _from_idx_to_inter;

   /// Weight derivative matrix
   std::vector<std::vector<double>> _der_matrix;
   /// Vector of indexes in the weight array corresponding to intervals boundaries (dim = num
   /// intervals + 1)
   std::vector<size_t> _delim_indexes;

   /// Vector of integrals of the weights
   vector_d _integral_weights;

   /// Number of weights
   size_t size = 0;
   /// Number of weights (as signed integer)
   i32 size_li = 0;
   /// Number of nodes
   size_t c_size = 0;
   /// Number of nodes (as signed integer)
   i32 c_size_li = 0;
};

template <typename O, typename T>
requires cpt::InterpolateCompatible<T, O>
O Discretize(const Grid1D &grid, const std::function<T(double)> &fnc,
             const std::function<O(size_t)> &make_zero_container)
{
   O res = make_zero_container(grid.c_size);
   for (size_t i = 0; i < grid.c_size; i++) {
      res[i] = fnc(grid._coord[i]);
   }
   return res;
}

} // namespace Interpolation