#ifndef INCLUDE_METHODS_INTERFACES_ADAPTIVE_SINGLE_STEP_METHOD_HPP_
#define INCLUDE_METHODS_INTERFACES_ADAPTIVE_SINGLE_STEP_METHOD_HPP_

#include <concepts>
#include "types.hpp"

namespace odelib {

/**
 * AdaptiveSingleStepMethod
 * An explicit, adaptive, single-step method for solving ODEs.
 */
template<typename Method>
concept AdaptiveSingleStepMethod = requires(Method met, double tol,
    double t, const Vectord<1>& x, double& h, const Vectord<1>& dv) {
  { Method::order() } -> std::same_as<int>;
  // { met.step(t, x, h, tol) } -> std::same_as<typeof(x), double>;
  { met.step(t, {}, h, tol) };
  // { met.step(t, x, h, tol, dv) } -> std::same_as<typeof(x), double>;
  // { met.hinted_step(t, {}, h, tol, {}) };
};

}  // namespace odelib

#endif  // INCLUDE_METHODS_INTERFACES_ADAPTIVE_SINGLE_STEP_METHOD_HPP_

