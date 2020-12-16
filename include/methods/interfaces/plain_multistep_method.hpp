#ifndef INCLUDE_METHODS_INTERFACES_PLAIN_MULTISTEP_METHOD_HPP_
#define INCLUDE_METHODS_INTERFACES_PLAIN_MULTISTEP_METHOD_HPP_

#include <concepts>
#include "problems/arenstorf.hpp"
#include "types.hpp"
#include "initial_value_problem.hpp"

namespace odelib {

/**
 * PlainMultistepMethod
 * An explicit, non-adaptive, multistep method for solving ODEs.
 */
template<typename Method, typename Dv = Arenstorf::Dv>
concept PlainMultistepMethod = requires(Method met, Dv f, double t,
    const Vectord<Dv::kDim>* x, double h, const Vectord<Dv::kDim>* dv) {
  { Method::order } -> std::same_as<const int&>;
  { Method::neededSteps } -> std::same_as<const int&>;
  { met.step(f, t, x, h, dv) } -> std::same_as<Vectord<Dv::kDim>>;
};

}  // namespace odelib

#endif  // INCLUDE_METHODS_INTERFACES_PLAIN_MULTISTEP_METHOD_HPP_

