#ifndef INCLUDE_METHODS_MOD_EULER_HPP_
#define INCLUDE_METHODS_MOD_EULER_HPP_

#include "initial_value_problem.hpp"
#include "types.hpp"

namespace odelib {

template<IvpDerivative D>
struct ModEuler {
  D f;

  static constexpr int order() {
    return 2;
  }

  inline Vectord<D::kDim> step(double t, const Vectord<D::kDim>& x,
                               double h) const {
    return hinted_step(t, x, h, f(t, x));
  }

  inline Vectord<D::kDim> hinted_step(double t, const Vectord<D::kDim>& x,
                              double h, const Vectord<D::kDim>& dv) const {
    return x + (dv*h + f(t+h, x+dv*h)*h)/2;
  }
};

}  // namespace odelib

#endif  // INCLUDE_METHODS_MOD_EULER_HPP_
