#ifndef INCLUDE_METHODS_INTERFACES_PLAIN_MULTISTEP_METHOD_HPP_
#define INCLUDE_METHODS_INTERFACES_PLAIN_MULTISTEP_METHOD_HPP_

#include <concepts>
#include "types.hpp"

namespace odelib {

/**
 * PlainMultistepMethod
 * An explicit, non-adaptive, multistep method for solving ODEs.
 */
template<typename Method>
concept PlainMultistepMethod = requires(Method met,
    double t, const Vectord<1>* x, double h, const Vectord<1>* dv) {
  { Method::order() } -> std::same_as<int>;
  { Method::neededSteps() } -> std::same_as<int>;
  // { met.step(t, x, h, dv) } -> std::same_as<Vectord<1>>;
  { met.step(t, nullptr, h, nullptr) };
};

}  // namespace odelib

#endif  // INCLUDE_METHODS_INTERFACES_PLAIN_MULTISTEP_METHOD_HPP_

