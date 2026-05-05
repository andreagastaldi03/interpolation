
#pragma once

#include <concepts>
#include "Interpolation/default.hh"

namespace Interpolation
{

using index_t = i32;

namespace cpt
{

// =============================================================================
// ===================== DISCRETIZATION CONCEPTS ===============================
// =============================================================================
/**
 * @brief Checks if the given pair of ReturnType and InnerType can be used in the interpolation
 * routine
 *
 * @tparam ReturnType The output type
 * @tparam InnerType The input type, typically an indexable container of ReturnType
 *
 * @note Albeit rare, InnerType needs not to be a container of ReturnType objects,
 * it suffice that InnerType is indexable and (InnerType[i] * w) can be summed to
 * the ReturnTpe object
 */
template <class ReturnType, class InnerType>
concept InterpolateCompatibleIndex
    = requires(ReturnType t, const InnerType &in, index_t i, double w) {
         { in[i] };
         { t += in[i] * w } -> std::same_as<ReturnType &>;
      } || requires(ReturnType t, const InnerType &in, size_t i, double w) {
         { in[i] };
         { t += in[i] * w } -> std::same_as<ReturnType &>;
      };

// "Element has .Evaluate(args...) and that result can be accumulated like before"
template <class ReturnType, class InnerType, class... Args>
concept InterpolateCompatibleEvaluate
    = requires(ReturnType t, const InnerType &in, index_t i, double w, Args &&...args) {
         { in[i].Evaluate(std::forward<Args>(args)...) };
         { t += in[i].Evaluate(std::forward<Args>(args)...) * w } -> std::same_as<ReturnType &>;
      } || requires(ReturnType t, const InnerType &in, size_t i, double w, Args &&...args) {
         { in[i].Evaluate(std::forward<Args>(args)...) };
         { t += in[i].Evaluate(std::forward<Args>(args)...) * w } -> std::same_as<ReturnType &>;
      };

template <class ReturnType, class InnerType, class... Args>
concept InterpolateCompatibleFunctor
    = requires(ReturnType t, const InnerType &in, index_t i, double w, Args &&...args) {
         { in[i](std::forward<Args>(args)...) };
         { t += in[i](std::forward<Args>(args)...) * w } -> std::same_as<ReturnType &>;
      } || requires(ReturnType t, const InnerType &in, std::size_t i, double w, Args &&...args) {
         { in[i](std::forward<Args>(args)...) };
         { t += in[i](std::forward<Args>(args)...) * w } -> std::same_as<ReturnType &>;
      };

template <class ReturnType, class InnerType, class... Args>
concept InterpolateCompatible = InterpolateCompatibleIndex<ReturnType, InnerType>
                             || InterpolateCompatibleEvaluate<ReturnType, InnerType, Args...>
                             || InterpolateCompatibleFunctor<ReturnType, InnerType, Args...>;

/**
 * @brief Checks if the container is indexable as a vector of vectors
 *
 * @tparam Container The container to be checked
 *
 * @warning Index types supported are only size_t and index_t
 */
template <class Container>
concept Index2DVecVec = requires(const Container &c, index_t i, index_t j) {
   { c[i][j] };
} || requires(const Container &c, index_t i, size_t j) {
   { c[i][j] };
} || requires(const Container &c, size_t i, index_t j) {
   { c[i][j] };
} || requires(const Container &c, size_t i, size_t j) {
   { c[i][j] };
};

/**
 * @brief Checks if the container is indexable as an Eigen::Matrix
 *
 * @tparam Container The container to be checked
 *
 * @warning Index types supported are only size_t and index_t
 */
template <class Container>
concept Index2DMat = requires(const Container &c, index_t i, index_t j) {
   { c(i, j) };
} || requires(const Container &c, index_t i, size_t j) {
   { c(i, j) };
} || requires(const Container &c, size_t i, index_t j) {
   { c(i, j) };
} || requires(const Container &c, size_t i, size_t j) {
   { c(i, j) };
};

template <class Container, class... Args>
concept Index2DVecVecEvaluate = requires(const Container &c, index_t i, index_t j, Args &&...args) {
   { c[i][j].Evaluate(std::forward<Args>(args)...) };
} || requires(const Container &c, index_t i, size_t j, Args &&...args) {
   { c[i][j].Evaluate(std::forward<Args>(args)...) };
} || requires(const Container &c, size_t i, index_t j, Args &&...args) {
   { c[i][j].Evaluate(std::forward<Args>(args)...) };
} || requires(const Container &c, size_t i, size_t j, Args &&...args) {
   { c[i][j].Evaluate(std::forward<Args>(args)...) };
};

template <class Container, class... Args>
concept Index2DMatEvaluate = requires(const Container &c, index_t i, index_t j, Args &&...args) {
   { c(i, j).Evaluate(std::forward<Args>(args)...) };
} || requires(const Container &c, index_t i, size_t j, Args &&...args) {
   { c(i, j).Evaluate(std::forward<Args>(args)...) };
} || requires(const Container &c, size_t i, index_t j, Args &&...args) {
   { c(i, j).Evaluate(std::forward<Args>(args)...) };
} || requires(const Container &c, size_t i, size_t j, Args &&...args) {
   { c(i, j).Evaluate(std::forward<Args>(args)...) };
};

template <class ReturnType, class InnerType>
concept InterpolateCompatible2DIndex
    = (Index2DMat<InnerType>  && (requires(ReturnType t, const InnerType &in, index_t i, index_t j, double w) {
         { t += in(i, j) * w } -> std::same_as<ReturnType &>;
      } || requires(ReturnType t, const InnerType &in, index_t i, size_t j, double w) {
         { t += in(i, j) * w } -> std::same_as<ReturnType &>;
      }|| requires(ReturnType t, const InnerType &in, size_t i, index_t j, double w) {
         { t += in(i, j) * w } -> std::same_as<ReturnType &>;
      }|| requires(ReturnType t, const InnerType &in, size_t i, size_t j, double w) {
         { t += in(i, j) * w } -> std::same_as<ReturnType &>;
      })) || (Index2DVecVec<InnerType>  && (requires(ReturnType t, const InnerType &in, index_t i, index_t j, double w) {
         { t += in[i][j] * w } -> std::same_as<ReturnType &>;
      } || requires(ReturnType t, const InnerType &in, index_t i, size_t j, double w) {
         { t += in[i][j] * w } -> std::same_as<ReturnType &>;
      }|| requires(ReturnType t, const InnerType &in, size_t i, index_t j, double w) {
         { t += in[i][j] * w } -> std::same_as<ReturnType &>;
      }|| requires(ReturnType t, const InnerType &in, size_t i, size_t j, double w) {
         { t += in[i][j] * w } -> std::same_as<ReturnType &>;
      })) ;

template <class ReturnType, class InnerType, class... Args>
concept InterpolateCompatible2DEvaluate
    = (Index2DMatEvaluate<InnerType, Args...>  && (requires(ReturnType t, const InnerType &in, index_t i, index_t j, double w, Args&&...args) {
         { t += in(i, j).Evaluate(std::forward<Args>(args)...) * w } -> std::same_as<ReturnType &>;
      } || requires(ReturnType t, const InnerType &in, index_t i, size_t j, double w, Args&&...args) {
         { t += in(i, j).Evaluate(std::forward<Args>(args)...) * w } -> std::same_as<ReturnType &>;
      }|| requires(ReturnType t, const InnerType &in, size_t i, index_t j, double w, Args&&...args) {
         { t += in(i, j).Evaluate(std::forward<Args>(args)...) * w } -> std::same_as<ReturnType &>;
      }|| requires(ReturnType t, const InnerType &in, size_t i, size_t j, double w, Args&&...args) {
         { t += in(i, j).Evaluate(std::forward<Args>(args)...) * w } -> std::same_as<ReturnType &>;
      })) || (Index2DVecVecEvaluate<InnerType, Args...>  && (requires(ReturnType t, const InnerType &in, index_t i, index_t j, double w, Args&&...args) {
         { t += in[i][j].Evaluate(std::forward<Args>(args)...) * w } -> std::same_as<ReturnType &>;
      } || requires(ReturnType t, const InnerType &in, index_t i, size_t j, double w, Args&&...args) {
         { t += in[i][j].Evaluate(std::forward<Args>(args)...) * w } -> std::same_as<ReturnType &>;
      }|| requires(ReturnType t, const InnerType &in, size_t i, index_t j, double w, Args&&...args) {
         { t += in[i][j].Evaluate(std::forward<Args>(args)...) * w } -> std::same_as<ReturnType &>;
      }|| requires(ReturnType t, const InnerType &in, size_t i, size_t j, double w, Args&&...args) {
         { t += in[i][j].Evaluate(std::forward<Args>(args)...) * w } -> std::same_as<ReturnType &>;
      })) ;

template <class ReturnType, class InnerType, class... Args>
concept InterpolateCompatible2D = InterpolateCompatible2DIndex<ReturnType, InnerType>
                               || InterpolateCompatible2DEvaluate<ReturnType, InnerType, Args...>;

/**
 * @brief Checks if the class K provides the 4 required static double->double functions:
 *        tis, tis_d, tps, tps_d
 *        that are required to be interpreted as a PIM struct (Physical Interpolation Map)
 *
 * @tparam PIM The type to be checked
 */
template <class K>
concept isPIM = requires(double x) {
   { K::tis(x) } -> std::same_as<double>; // to interpolate space
   { K::tis_d(x) } -> std::same_as<double>;
   { K::tps(x) } -> std::same_as<double>; // to physical space
   { K::tps_d(x) } -> std::same_as<double>;
};

} // namespace cpt
} // namespace Interpolation
