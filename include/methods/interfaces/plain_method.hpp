#ifndef INCLUDE_METHODS_INTERFACES_PLAIN_METHOD_HPP_
#define INCLUDE_METHODS_INTERFACES_PLAIN_METHOD_HPP_

#include <concepts>
#include "types.hpp"

namespace odelib {

/**
 * PlainMethod
 * An explicit, non-adaptive, single-step method for solving ODEs.
 */
template<typename Method>
concept PlainMethod = std::default_initializable<Method> && requires(Method met,
    double t, const Vectord<1>& x, double h, const Vectord<1>& dv) {
  { Method::order() } -> std::same_as<int>;
  // { met.step(t, x, h) } -> std::same_as<typeof(x)>;
  { met.step(t, {}, h) };
  // { met.step(t, x, h, dv) } -> std::same_as<typeof(x)>;
  { met.hinted_step(t, {}, h, {}) };
};

}  // namespace odelib

#endif  // INCLUDE_METHODS_INTERFACES_PLAIN_METHOD_HPP_

