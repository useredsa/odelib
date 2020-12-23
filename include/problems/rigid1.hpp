#ifndef INCLUDE_PROBLEMS_RIGID_1D_HPP_
#define INCLUDE_PROBLEMS_RIGID_1D_HPP_

#include <cmath>
#include "types.hpp"

namespace odelib {

/**
 * A rigid/stiff problem of dimension 1.
 * Contains the partial derivative with respect to space of the function,
 * allowing to solve the problem using Newton's method.
 */
struct Rigid1 {
  static inline double t0() { return 0; }
  static inline Vectord<1> x0() { return Vectord<1>{-1}; }

  struct Dv {
    static constexpr int kDim = 1;

    inline Vectord<1> operator()(double t, const Vectord<1>& x) const {
      return Vectord<1>{5*exp(5*t)*(x[0] - t)*(x[0] - t) + 1};
    }

    inline Vectord<1> pdvx(double t, const Vectord<1>& x) const {
      return Vectord<1>{10*exp(5*t)*(x[0] - t)};
    }
  };

  static inline Vectord<1> analyticalSol(double t) {
    return Vectord<1>{t - exp(-5 * t)};
  }
};

}  // namespace odelib

#endif  // INCLUDE_PROBLEMS_RIGID_1D_HPP_
