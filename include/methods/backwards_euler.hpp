#ifndef INCLUDE_METHODS_BACKWARDS_EULER_HPP_
#define INCLUDE_METHODS_BACKWARDS_EULER_HPP_

#include "methods/interfaces/implicit_equation.hpp"
#include "initial_value_problem.hpp"
#include "types.hpp"

namespace odelib {

/**
 * Backward Euler's method.
 * The simplest implicit method.
 * Gives an approximation of order 1.
 * 
 * It can also be seen as the backward differentiation formula of order 1.
 * 
 * x_{n+1} = x_n + h*(t_{n+1}, x_{n+1})
 */
struct BackwardsEuler {
  static constexpr int kOrder = 1;

// Hinted equation is equal to actual equation
// because the derivative is not used

  template <IvpDerivative D>
  inline LfImplicitEquation<D> equation(const D& f, double t,
      const Vectord<D::kDim>& x, double h) const {
    return LfImplicitEquation<D>(x, h, t+h);
  }

  template <IvpDerivative D>
  inline LfImplicitEquation<D> hinted_equation(const D& f, double t,
      const Vectord<D::kDim>& x, double h, const Vectord<D::kDim>& dv) const {
    return LfImplicitEquation<D>(x, h, t+h);
  }
};

}  // namespace odelib

#endif  // INCLUDE_METHODS_BACKWARDS_EULER_HPP_
