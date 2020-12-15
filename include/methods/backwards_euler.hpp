#ifndef INCLUDE_METHODS_BACKWARDS_EULER_HPP_
#define INCLUDE_METHODS_BACKWARDS_EULER_HPP_

#include "implicit_equation.hpp"
#include "initial_value_problem.hpp"
#include "types.hpp"

namespace odelib {

// Hinted equation is equal to actual equation
// because the derivative is not used

template<IvpDerivative D>
struct BackwardsEuler {
  D f;

  static constexpr int order() {
    return 1;
  }

  inline LfImplicitEquation<D> equation(double t,
                     const Vectord<D::kDim>& x, double h) {
    return LfImplicitEquation<D>(x, h, t+h);
  }

  inline LfImplicitEquation<D> hinted_equation(double t,
      const Vectord<D::kDim>& x, double h, const Vectord<D::kDim>& dv) {
    return LfImplicitEquation<D>(x, h, t+h);
  }
};

}  // namespace odelib

#endif  // INCLUDE_METHODS_BACKWARDS_EULER_HPP_
