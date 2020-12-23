#ifndef INCLUDE_METHODS_FEHLBERG_HPP_
#define INCLUDE_METHODS_FEHLBERG_HPP_

#include <algorithm>
#include <utility>
#include "types.hpp"

namespace odelib {

/**
 * Fehlberg's Method
 * An adaptive Runge-Kutta method
 * 
 * Uses two Runge-Kutta methods of order's 4 and 5
 * which share intermeditate steps.
 * Thanks to this, with only 6 evaluations we can compute
 * a fifth order Runge-Kutta method and
 * an estimate of the error.
 */
struct Fehlberg {
  static constexpr int kOrder = 4;

  template <IvpDerivative D>
  inline std::pair<Vectord<D::kDim>, double> step(D f, double t,
      const Vectord<D::kDim>& x, double& h, double tolerance) const {
    Vectord<D::kDim> k[6];
    k[0] = h*f(t, x);
    k[1] = h*f(t +   1/4.0*h, x + 1/4.0*k[0]);
    k[2] = h*f(t +   3/8.0*h, x + (3*k[0] + 9*k[1])/32);
    k[3] = h*f(t + 12/13.0*h, x + (1932*k[0] - 7200*k[1] + 7296*k[2])/2197);
    k[4] = h*f(t +         h, x + 439/216.0*k[0] - 8*k[1] + 3680/513.0*k[2]
                                - 845/4104.0*k[3]);
    k[5] = h*f(t +   1/2.0*h, x - 8/27.0*k[0] + 2*k[1] - 3544/2565.0*k[2]
                                + 1859/4104.0*k[3] - 11/40.0*k[4]);

    Vectord<D::kDim> rk4 = x + 25/216.0*k[0] + 1408/2565.0*k[2]
                          + 2197/4104.0*k[3] - 1/5.0*k[4];
    Vectord<D::kDim> rk5 = x + 16/135.0*k[0] + 6656/12825.0*k[2]
                        + 28561/56430.0*k[3] - 9/50.0*k[4] + 2/55.0*k[5];
    double error = (rk5 - rk4).norm();
    double q = std::pow(tolerance * h / (2*error), 0.25);
    q = std::max(0.1, std::min(q, 4.0));
    h *= q;
    return {rk5, error};
  }
};

}  // namespace odelib

#endif  // INCLUDE_METHODS_FEHLBERG_HPP_
