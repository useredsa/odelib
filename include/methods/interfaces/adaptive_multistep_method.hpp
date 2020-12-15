#ifndef INCLUDE_METHODS_INTERFACES_ADAPTIVE_MULTISTEP_METHOD_HPP_
#define INCLUDE_METHODS_INTERFACES_ADAPTIVE_MULTISTEP_METHOD_HPP_

#include <concepts>
#include "types.hpp"

namespace odelib {

/**
 * AdaptiveMultistepMethod
 * An explicit, adaptive, multistep method for solving ODEs.
 */
template<typename Method>
concept AdaptiveMultistepMethod = requires(Method met, double tol,
    double t, const Vectord<1>* x, double& h, const Vectord<1>* dv) {
  { Method::order() } -> std::same_as<int>;
  { Method::neededSteps() } -> std::same_as<int>;
  // { met.step(t, x, h, dv, tol) }
  //     -> std::same_as<std::pair<Vectord<1>, double>>;
  { met.step(t, nullptr, h, nullptr, tol) };
};

}  // namespace odelib

#endif  // INCLUDE_METHODS_INTERFACES_ADAPTIVE_MULTISTEP_METHOD_HPP_

