#ifndef INCLUDE_PROBLEMS_RIGID_1D_HPP_
#define INCLUDE_PROBLEMS_RIGID_1D_HPP_

#include <cmath>
#include "Eigen/Dense"
// #include "initial_value_problem.hpp"

using Eigen::Vectord;

constexpr int kDim = 1;

// class Rigid_1d : InitialValueProblem<kDim> {
//   Rigid_1d(double t, const Vectord<N> &x)
//     : InitialValueProblem(t, x) {}
//
//   inline Vectord<kDim> derivative(double t, const Vectord<kDim>& x) {
//     return Vectord<1>{
//       5 * exp(5*t) * (x[0] - t) * (x[0] - t) + 1
//     };
//   }
// }

struct Rigid_1d {
  double t0;
  Vectord<kDim> x0;

  Rigid_1d(double t0, Vectord<kDim> x0)
    : t0(t0), x0(x0) {}

  inline double get_t0() {
    return t0;
  }

  inline Vectord<kDim> get_x0() {
    return x0;
  }

  inline Vectord<kDim> derivative(double t, const Vectord<kDim>& x) {
    return Vectord<1>{
      5 * exp(5*t) * (x[0] - t) * (x[0] - t) + 1
    };
  }

  inline Vectord<kDim> dy_derivative(double t, const Vectord<kDim>& x) {
    return Vectord<1>{
      10 * exp(5*t) * (x[0] - t)
    };
  }
};

struct analytical_sol {
  inline Vectord<kDim> operator()(double t) {
    return Vectord<1>{
      t - exp(-5 * t)
    };
  }
};


#endif  // INCLUDE_PROBLEMS_RIGID_1D_HPP_
