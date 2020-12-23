#ifndef INCLUDE_METHODS_INTERFACES_PLAIN_ADAPTIVE_METHOD_HPP_
#define INCLUDE_METHODS_INTERFACES_PLAIN_ADAPTIVE_METHOD_HPP_

#include <concepts>
#include <utility>
#include "problems/arenstorf.hpp"
#include "types.hpp"

namespace odelib {

/**
 * PlainAdaptiveMethod
 * An explicit, adaptive, single-step method for solving ODEs.
 */
template <typename Method, typename Dv = Arenstorf::Dv>
concept PlainAdaptiveMethod = requires(Method met, Dv f, double t,
    const Vectord<Dv::kDim>& x, double& h, const Vectord<Dv::kDim>& dv,
    double tol) {
  { Method::kOrder } -> std::same_as<const int&>;
  { met.step(f, t, x, h, tol) }
      -> std::same_as<std::pair<Vectord<Dv::kDim>, double>>;
  // { met.hinted_step(f, t, x, h, tol, dv) }
  //     -> std::same_as<Vectord<Dv::kDim>, double>;
};

}  // namespace odelib

#endif  // INCLUDE_METHODS_INTERFACES_PLAIN_ADAPTIVE_METHOD_HPP_
