#ifndef INCLUDE_METHODS_SSPRK3_HPP_
#define INCLUDE_METHODS_SSPRK3_HPP_

#include "initial_value_problem.hpp"
#include "types.hpp"

namespace odelib {

/**
 * Third Order Strong Stability Preserving Explicit Runge-Kutta
 * 
 * Tableau:
 * 
 *  0  ┃
 *  1  ┃  1
 * 1/2 ┃ 1/4 1/4
 * ━━━━╋━━━━━━━━━━━━━
 *     ┃ 1/6 1/6 2/3
 */
struct SSPRK3 {
  static constexpr int kOrder = 3;

  template <IvpDerivative D>
  inline Vectord<D::kDim> step(D f, double t, const Vectord<D::kDim>& x,
      double h) const {
    return hinted_step(f, t, x, h, f(t, x));
  }

  template <IvpDerivative D>
  inline Vectord<D::kDim> hinted_step(D f, double t, const Vectord<D::kDim>& x,
      double h, const Vectord<D::kDim>& dv) const {
    Vectord<D::kDim> k[3];
    k[0] = h*dv;
    k[1] = h*f(t + h, x + k[0]);
    k[2] = h*f(t + h/2, x + k[0]/4 + k[1]/4);
    return x + (k[0] + k[1] + k[2]*4)/6;
  }
};

}  // namespace odelib

#endif  // INCLUDE_METHODS_SSPRK3_HPP_

