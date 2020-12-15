#ifndef INCLUDE_METHODS_RK4_HPP_
#define INCLUDE_METHODS_RK4_HPP_

#include "initial_value_problem.hpp"
#include "types.hpp"

namespace odelib {

template<IvpDerivative D>
struct RK4 {
  D f;

  static constexpr int order() {
    return 4;
  }

  inline Vectord<D::kDim> step(double t, const Vectord<D::kDim>& x,
                               double h) const {
    return hinted_step(t, x, h, f(t, x));
  }

  inline Vectord<D::kDim> hinted_step(double t, const Vectord<D::kDim>& x,
                               double h, const Vectord<D::kDim>& d) const {
    Vectord<D::kDim> k[4];
    k[0] = d*h;
    k[1] = f(t+h/2, x + k[0]/2)*h;
    k[2] = f(t+h/2, x + k[1]/2)*h;
    k[3] = f(t+h, x + k[2])*h;
    return (k[0] + 2*k[1] + 2*k[2] + k[3])/6;
  }
};

}  // namespace odelib

#endif  // INCLUDE_METHODS_RK4_HPP_
