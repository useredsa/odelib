#ifndef INCLUDE_METHODS_TRAPEZOIDAL_HPP_
#define INCLUDE_METHODS_TRAPEZOIDAL_HPP_

#include "implicit_equation.hpp"
#include "initial_value_problem.hpp"
#include "types.hpp"

namespace odelib {

struct Trapezoidal {
  static constexpr int order = 2;

  template<IvpDerivative D>
  inline LfImplicitEquation<D> equation(D f, double t,
      const Vectord<D::kDim>& x, double h) {
    return hinted_equation(t, x, h, f(t, x));
  }

  template<IvpDerivative D>
  inline LfImplicitEquation<D> hinted_equation(D f, double t,
      const Vectord<D::kDim>& x, double h, const Vectord<D::kDim>& dv) {
    return LfImplicitEquation<D>(x + h/2*dv, h/2, t+h);
  }
};

}  // namespace odelib

#endif  // INCLUDE_METHODS_TRAPEZOIDAL_HPP_
