#ifndef INCLUDE_METHODS_EULER_HPP_
#define INCLUDE_METHODS_EULER_HPP_

#include "initial_value_problem.hpp"
#include "types.hpp"

namespace odelib {

struct Euler {
  static constexpr int order = 1;

  template<IvpDerivative D>
  inline Vectord<D::kDim> step(D f, double t, const Vectord<D::kDim>& x,
      double h) const {
    return hinted_step(f, t, x, h, f(t, x));
  }

  template<IvpDerivative D>
  inline Vectord<D::kDim> hinted_step(D f, double t, const Vectord<D::kDim>& x,
      double h, const Vectord<D::kDim>& dv) const {
    return x + dv*h;
  }
};

}  // namespace odelib

#endif  // INCLUDE_METHODS_EULER_HPP_
