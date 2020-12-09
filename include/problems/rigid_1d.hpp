#ifndef INCLUDE_PROBLEMS_RIGID_1D_HPP_
#define INCLUDE_PROBLEMS_RIGID_1D_HPP_

#include <cmath>
#include "Eigen/Dense"

namespace rigid_1d {

using Eigen::Vectord;
using std::sqrt;

constexpr int kDim = 1;

struct derivative {

  inline Vectord<kDim> operator()(double t, const Vectord<kDim>& x) {
    return Vectord<1>{
      5 * exp(5*t) * (x[0] - t) * (x[0] - t) + 1
    };
  }

  inline Vectord<kDim> d_y(double t, const Vectord<kDim>& x) {
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

}  // namespace rigid_1d

#endif  // INCLUDE_PROBLEMS_RIGID_1D_HPP_
