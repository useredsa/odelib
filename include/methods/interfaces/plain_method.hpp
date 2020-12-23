#ifndef INCLUDE_METHODS_INTERFACES_PLAIN_METHOD_HPP_
#define INCLUDE_METHODS_INTERFACES_PLAIN_METHOD_HPP_

#include <concepts>
#include "initial_value_problem.hpp"
#include "problems/taylor1.hpp"
#include "types.hpp"

namespace odelib {

/**
 * PlainMethod
 * An explicit, non-adaptive, single-step method for solving ODEs.
 */
template <typename Method, typename Dv = Taylor1::Dv>
concept PlainMethod = std::default_initializable<Method> && requires(Method met,
    Dv f, double t, const Vectord<Dv::kDim>& x, double h,
    const Vectord<Dv::kDim>& dv) {
  { Method::kOrder } -> std::same_as<const int&>;
  { met.step(f, t, x, h) } -> std::same_as<Vectord<Dv::kDim>>;
  { met.hinted_step(f, t, x, h, dv) } -> std::same_as<Vectord<Dv::kDim>>;
};

}  // namespace odelib

#endif  // INCLUDE_METHODS_INTERFACES_PLAIN_METHOD_HPP_

