#ifndef INCLUDE_METHODS_INTERFACES_BACKWARD_DIFFERENTIATION_FORMULA_HPP_
#define INCLUDE_METHODS_INTERFACES_BACKWARD_DIFFERENTIATION_FORMULA_HPP_

#include "types.hpp"
#include "implicit_equation.hpp"

namespace odelib {

/**
 * BackwardDifferentiationFormula
 * Special family of implicit, non-adaptive, multistep methods.
 */
template<typename Method>
concept BackwardDifferentiationFormula = requires(Method met,
    double t, const Vectord<1>* x, double h, const Vectord<1>& dv) {
  { Method::order() } -> std::same_as<int>;
  { Method::neededSteps() } -> std::same_as<int>;
  // { ImplicitFormula<decltype(met.step(t, x, h))> };
  // { ImplicitEquation<decltype(met.step(t, nullptr, h))> };
};

}  // namespace odelib

#endif // INCLUDE_METHODS_INTERFACES_BACKWARD_DIFFERENTIATION_FORMULA_HPP_

