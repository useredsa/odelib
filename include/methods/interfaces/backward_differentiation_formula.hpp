#ifndef INCLUDE_METHODS_INTERFACES_BACKWARD_DIFFERENTIATION_FORMULA_HPP_
#define INCLUDE_METHODS_INTERFACES_BACKWARD_DIFFERENTIATION_FORMULA_HPP_

#include <concepts>
#include "problems/arenstorf.hpp"
#include "types.hpp"
#include "implicit_equation.hpp"

namespace odelib {

/**
 * BackwardDifferentiationFormula
 * Special family of implicit, non-adaptive, multistep methods.
 */
template <typename Method, typename Dv = Arenstorf::Dv>
concept BackwardDifferentiationFormula = requires(Method met, Dv f, double t,
    const Vectord<Dv::kDim>* x, double h, const Vectord<Dv::kDim>& dv) {
  { Method::kOrder } -> std::same_as<const int&>;
  { Method::kNeededSteps } -> std::same_as<const int&>;
  { met.equation(f, t, x, h) } -> ImplicitEquation;
};

}  // namespace odelib

#endif  // INCLUDE_METHODS_INTERFACES_BACKWARD_DIFFERENTIATION_FORMULA_HPP_
