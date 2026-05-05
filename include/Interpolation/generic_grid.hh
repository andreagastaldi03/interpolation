#pragma once

#include "Interpolation/default.hh"

namespace Interpolation
{

namespace Generic
{

/**
 * @brief Generic standardized one-dimensional grid
 *
 * Stores the points, the weights, the derivative matrix and the methods to interpolate
 */
struct StandardGrid {

   enum class STRATEGY : unsigned {
      NAIVE = 0,
      FBF   = 1,
      SBF   = 2,
   };

   /**
    * @brief Construct a new Standard Grid object from an input vector
    *
    * @param input The vector of nodes to be used
    *
    * @note The input vector is checked and sorted. If extremes points are not
    * [-1, 1], they are appended.
    */
   StandardGrid(const vector_d &input);

   /**
    * @brief Construct a new Standard Grid object from an input routine
    *
    * @param fnc The routine to produce the grid points, as function of the
    *            index and the total number of points
    * @param N The degree of the interpolating polynomial (total # of points = p+1)
    *
    * @note It expects that the lower bound is -1 and the upper bound +1
    */
   StandardGrid(const std::function<double(size_t, size_t)> &fnc, size_t p);

   /**
    * @brief Interpolate a view of a vector on the grid
    *
    * @param t        The point in which we want the interpolation
    * @param fj       The vector of values to interpolate
    * @param start    The begin of the window in the vector.
    * @param end      The end of the window in the vector (inclusive)
    * @param str      The interpolation strategy to be used (default: first barycentric formula)
    * @return         The interpolated value
    */
   double interpolate(double t, const vector_d &fj, size_t start, size_t end,
                      STRATEGY str = STRATEGY::FBF) const;
   /**
    * @brief Interpolate the derivative of a view of a vector on the grid
    *
    * @param t        The point in which we want the interpolation
    * @param fj       The vector of values to interpolate
    * @param start    The begin of the window in the vector.
    * @param end      The end of the window in the vector (inclusive)
    * @param str      The interpolation strategy to be used (default: first barycentric formula)
    * @return         The interpolated value
    */
   double interpolate_der(double t, const vector_d &fj, size_t start, size_t end,
                          STRATEGY str = STRATEGY::FBF) const;

   /**
    * @brief Interpolate the derivative of a view of a vector on the grid
    *
    * @param t        The point in which we want the interpolation
    * @param fj       The vector of values to interpolate
    * @param start    The begin of the window in the vector.
    * @param end      The end of the window in the vector (inclusive)
    * @param str      The interpolation strategy to be used (default: first barycentric formula)
    * @return         The interpolated value
    *
    * @note This version takes the input vector, applies the derivative matrix,
    * and then calls the standard interpolation.
    */
   double interpolate_der_v2(double t, const vector_d &fj, size_t start, size_t end,
                             STRATEGY str = STRATEGY::FBF) const;

   /**
    * @brief Polynomial weight at index j, naive expression
    *
    * @param t The point in which to evaluate the weight
    * @param j The index of the weight to evaluate
    * @return  The value of the weight at the point
    */
   double poli_weight(double t, size_t j) const;

   /**
    * @brief Polynomial weight at index j for first barycentric formula
    *
    * @param t The point in which to evaluate the weight
    * @param j The index of the weight to evaluate
    * @return  The value of the weight at the point
    */
   double poli_weight_fbf(double t, size_t j) const;

   /**
    * @brief Polynomial weight at index j for first barycentric formula
    *
    * @param t     The point in which to evaluate the weight
    * @param j     The index of the weight to evaluate
    * @param monic The value of the monic Lagrange polynomial \f$ l(x) \f$
    * @return  The value of the weight at the point
    */
   double poli_weight_fbf(double t, size_t j, double monic) const;

   /**
    * @brief Polynomial weight at index j for first barycentric formula
    *
    * @param t The point in which to evaluate the weight
    * @param j The index of the weight to evaluate
    * @return  The value of the weight at the point
    */
   double poli_weight_sbf(double t, size_t j) const;

   /**
    * @brief Polynomial weight at index j for first barycentric formula
    *
    * @param t     The point in which to evaluate the weight
    * @param j     The index of the weight to evaluate
    * @param monic The denominator for the second barycentric formula
    * @return  The value of the weight at the point
    */
   double poli_weight_sbf(double t, size_t j, double den) const;

   /**
    * @brief Derivative of the polynomial weight at index j
    *
    * @param t The point in which to evaluate the derivative weight
    * @param j The index of the weight to evaluate
    * @return  The value of the derivative of the weight at the point (i.e. \f$ \partial w_j(t)
    * /
    * \partial t\f$)
    */
   double poli_weight_der(double t, size_t j) const;
   /// Same, for first barycentric formula
   double poli_weight_fbf_der(double t, size_t j) const;
   /// Same, for first barycentric formula with precomputed monic poly
   double poli_weight_fbf_der(double t, size_t j, double monic) const;
   /// Same, for second barycentric formula
   double poli_weight_sbf_der(double t, size_t j) const;
   /// Same, for second barycentric formula with precomputed denominator
   double poli_weight_sbf_der(double t, size_t j, double den) const;

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

   /// Degree of interpolating polynomial (# of points = p+1)
   size_t _p;
   /// Grid nodes
   vector_d _tj;
   /// \f$ \lambda_j \f$ weights for barycentric formulas
   vector_d _lambdaj;
   /// The derivative matrix
   std::vector<vector_d> _Dij;
};

} // namespace Generic
} // namespace Interpolation
