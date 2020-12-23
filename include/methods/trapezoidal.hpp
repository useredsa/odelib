#ifndef INCLUDE_METHODS_TRAPEZOIDAL_HPP_
#define INCLUDE_METHODS_TRAPEZOIDAL_HPP_

#include "methods/interfaces/implicit_equation.hpp"
#include "initial_value_problem.hpp"
#include "types.hpp"

namespace odelib {

/**
 * Trapezoidal Method or Trapezoidal Rule.
 * A fixed step implicit method of order 2.
 * 
 * x_{n+1} = x_n + h/2*(f(t_n, x_n) + f(t_{n+1}, x_{n+1}))
 */
struct Trapezoidal {
  static constexpr int kOrder = 2;

  template <IvpDerivative D>
  inline LfImplicitEquation<D> equation(const D& f, double t,
      const Vectord<D::kDim>& x, double h) const {
    return hinted_equation(f, t, x, h, f(t, x));
  }

  template <IvpDerivative D>
  inline LfImplicitEquation<D> hinted_equation(const D& f, double t,
      const Vectord<D::kDim>& x, double h, const Vectord<D::kDim>& dv) const {
    return LfImplicitEquation<D>(x + h/2*dv, h/2, t+h);
  }
};

}  // namespace odelib

#endif  // INCLUDE_METHODS_TRAPEZOIDAL_HPP_
