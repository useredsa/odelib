#ifndef INCLUDE_METHODS_INTERFACES_PLAIN_IMPLICIT_HPP_
#define INCLUDE_METHODS_INTERFACES_PLAIN_IMPLICIT_HPP_

#include <concepts>
#include "problems/arenstorf.hpp"
#include "types.hpp"
#include "implicit_equation.hpp"

namespace odelib {

/**
 * PlainImplicitMethod
 * An implicit, non-adaptive, single-step method for solving ODEs.
 */
template <typename Method, typename Dv = Arenstorf::Dv>
concept PlainImplicitMethod = requires(Method met, Dv f, double t,
    const Vectord<Dv::kDim>& x, double h, const Vectord<Dv::kDim>& dv) {
  { Method::kOrder } -> std::same_as<const int&>;
  { met.equation(f, t, x, h) } -> ImplicitEquation;
  { met.hinted_equation(f, t, x, h, dv) } -> ImplicitEquation;
};

}  // namespace odelib

#endif  // INCLUDE_METHODS_INTERFACES_PLAIN_IMPLICIT_HPP_
