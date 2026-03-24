#pragma once

#include "Interpolation/default.hh"

namespace Interpolation
{
namespace Chebyshev
{

/**
 * @brief Standardized one-dimensional Chebyshev grid, see
 * https://doi.org/10.48550/arXiv.2112.09703 https://doi.org/10.1140/epjc/s10052-022-10223-1
 *
 * Stores the points, the weights, the derivative matrix and the methods to interpolate
 */
struct StandardGrid {
   /**
    * @brief Construct a new Standard Grid object
    *
    * @param p The degree of the interpolating polynomial
    * (the number of points of the grid is p+1).
    */
   StandardGrid(size_t p = 3);

   /**
    * @brief Interpolate a view of a vector on the grid
    *
    * @param t        The point in which we want the interpolation
    * @param fj       The vector of values to interpolate
    * @param start    The begin of the window in the vector.
    * @param end      The end of the window in the vector (inclusive)
    * @return         The interpolated value
    */
   double interpolate(double t, const vector_d &fj, size_t start, size_t end) const;

   /**
    * @brief Interpolate the derivative of a view of a vector on the grid
    *
    * @param t        The point in which we want the interpolation
    * @param fj       The vector of values to interpolate
    * @param start    The begin of the window in the vector.
    * @param end      The end of the window in the vector (inclusive)
    * @return         The interpolated value
    */
   double interpolate_der(double t, const vector_d &fj, size_t start, size_t end) const;

   /**
    * @brief Polynomial weight at index j
    *
    * @param t The point in which to evaluate the weight
    * @param j The index of the weight to evaluate
    * @return  The value of the weight at the point
    */
   double poli_weight(double t, size_t j) const;

   /**
    * @brief Polynomial weight at index j
    *
    * @param t   The point in which to evaluate the weight
    * @param j   The index of the weight to evaluate
    * @param den The denominator for the second barycentric formula
    * @return    The value of the weight at the point
    *
    * This version is used when @interpolate is called, whereas previous
    * version is useful for other applications
    */
   double poli_weight(double t, size_t j, double den) const;

   /**
    * @brief Derivative of the polynomial weight at index j
    *
    * @param t The point in which to evaluate the derivative weight
    * @param j The index of the weight to evaluate
    * @return  The value of the derivaive of the weight at the point (i.e. \f$ \partial w_j(t)
    * /
    * \partial t\f$)
    */
   double poli_weight_der(double t, size_t j) const;

   /**
    * @brief Derivative of the polynomial weight at index j
    *
    * @param t   The point in which to evaluate the derivative weight
    * @param j   The index of the weight to evaluate
    * @param den The denominator for the interpolan
    * @return    The value of the derivaive of the weight at the point (i.e. \f$ \partial w_j(t)
    * /
    * \partial t\f$)
    */
   double poli_weight_der(double t, size_t j, double den) const;

   /**
    * @brief Applies the derivative matrix to the input in-place
    *
    * @param fj    Input view of discretized values, it is modified!
    * @param start The begin of the window in the vector.
    * @param end   The end of the window in the vector (inclusive)
    */
   void apply_D(vector_d &fj, size_t start, size_t end) const;

   /**
    * @brief Discretize the given function on the grid
    *
    * @param fnc       The function to be discretized
    * @return vector_d The vector of discretized values
    */
   vector_d discretize(const std::function<double(double)> &fnc) const;

   /**
    * @brief    Get the j-th grid point
    *
    * @param j  The index on the grid
    * @return   The value of the grid at the asked index
    */
   double t(size_t j) const
   {
      return _tj[j];
   }

   /// Degree of interpolating polynomial (# of points = N+1)
   size_t _p;
   /// Grid nodes
   vector_d _tj;
   /// \f$ \beta_j \f$ values, see right below Eq. 2.5 of
   /// https://doi.org/10.48550/arXiv.2112.09703
   vector_d _betaj;
   /// The derivative matrix
   std::vector<vector_d> _Dij;
};

} // namespace Chebyshev
} // namespace Interpolation
