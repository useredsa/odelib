#ifndef INCLUDE_METHODS_BACKWARDS_EULER_HPP_
#define INCLUDE_METHODS_BACKWARDS_EULER_HPP_

#include "implicit_equation.hpp"
#include "initial_value_problem.hpp"
#include "types.hpp"

namespace odelib {

// Hinted equation is equal to actual equation
// because the derivative is not used

struct BackwardsEuler {
  static constexpr int order = 1;

  template<IvpDerivative D>
  inline LfImplicitEquation<D> equation(D f, double t,
      const Vectord<D::kDim>& x, double h) {
    return LfImplicitEquation<D>(f, x, h, t+h);
  }

  template<IvpDerivative D>
  inline LfImplicitEquation<D> hinted_equation(D f, double t,
      const Vectord<D::kDim>& x, double h, const Vectord<D::kDim>& dv) {
    return LfImplicitEquation<D>(x, h, t+h);
  }
};

}  // namespace odelib

#endif  // INCLUDE_METHODS_BACKWARDS_EULER_HPP_
