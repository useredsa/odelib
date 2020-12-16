#ifndef INCLUDE_METHODS_MOD_EULER_HPP_
#define INCLUDE_METHODS_MOD_EULER_HPP_

#include "initial_value_problem.hpp"
#include "types.hpp"

namespace odelib {

struct ModEuler {
  static constexpr int order = 2;

  template<IvpDerivative D>
  inline Vectord<D::kDim> step(D f, double t, const Vectord<D::kDim>& x,
      double h) const {
    return hinted_step(t, x, h, f(t, x));
  }

  template<IvpDerivative D>
  inline Vectord<D::kDim> hinted_step(D f, double t, const Vectord<D::kDim>& x,
      double h, const Vectord<D::kDim>& dv) const {
    return x + (dv*h + f(t+h, x+dv*h)*h)/2;
  }
};

}  // namespace odelib

#endif  // INCLUDE_METHODS_MOD_EULER_HPP_
