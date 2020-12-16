#ifndef INCLUDE_METHODS_INTERFACES_ADAPTIVE_MULTISTEP_METHOD_HPP_
#define INCLUDE_METHODS_INTERFACES_ADAPTIVE_MULTISTEP_METHOD_HPP_

#include <concepts>
#include "problems/arenstorf.hpp"
#include "types.hpp"

namespace odelib {

/**
 * AdaptiveMultistepMethod
 * An explicit, adaptive, multistep method for solving ODEs.
 */
template<typename Method, typename Dv = Arenstorf::Dv>
concept AdaptiveMultistepMethod = requires(Method met, Dv f, double t,
    const Vectord<Dv::kDim>* x, double& h, const Vectord<Dv::kDim>* dv,
    double tol) {
  { Method::order } -> std::same_as<const int&>;
  { Method::neededSteps } -> std::same_as<const int&>;
  { met.step(f, t, x, h, dv, tol) }
      -> std::same_as<std::pair<Vectord<Dv::kDim>, double>>;
};

}  // namespace odelib

#endif  // INCLUDE_METHODS_INTERFACES_ADAPTIVE_MULTISTEP_METHOD_HPP_

