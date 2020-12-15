#ifndef INCLUDE_METHODS_INTERFACES_PLAIN_IMPLICIT_HPP_
#define INCLUDE_METHODS_INTERFACES_PLAIN_IMPLICIT_HPP_

#include "types.hpp"
#include "implicit_equation.hpp"

namespace odelib {

/**
 * PlainImplicitMethod
 * An implicit, non-adaptive, single-step method for solving ODEs.
 */
template<typename Method>
concept PlainImplicitMethod = requires(Method met,
    double t, const Vectord<1>& x, double h, const Vectord<1>& dv) {
  { Method::order() } -> std::same_as<int>;
  // { ImplicitFormula<decltype(met.equation(t, x, h))> };
  { ImplicitEquation<decltype(met.equation(t, {}, h))> };
  // { ImplicitFormula<decltype(met.hinted_equation(t, x, h, dv))> };
  { ImplicitEquation<decltype(met.hinted_equation(t, {}, h, dv))> };
};

}  // namespace odelib

#endif // INCLUDE_METHODS_INTERFACES_PLAIN_IMPLICIT_HPP_

