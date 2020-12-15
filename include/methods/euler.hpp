#ifndef INCLUDE_METHODS_EULER_HPP_
#define INCLUDE_METHODS_EULER_HPP_

#include "initial_value_problem.hpp"
#include "types.hpp"

namespace odelib {

template<IvpDerivative D>
struct Euler {
  D f;

  static constexpr int order() {
    return 1;
  }

  inline Vectord<D::kDim> step(double t, const Vectord<D::kDim>& x,
                               double h) const {
    return hinted_step(t, x, h, f(t, x));
  }

  inline Vectord<D::kDim> hinted_step(double t, const Vectord<D::kDim>& x,
                              double h, const Vectord<D::kDim>& dv) const {
    return x + dv*h;
  }
};

}  // namespace odelib

#endif  // INCLUDE_METHODS_EULER_HPP_
