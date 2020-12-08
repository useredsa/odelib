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
  };

  struct partial2 {
    inline Vectord<kDim> operator()(double t, const Vectord<kDim>& x) {
      double e = exp(5*t);
      double diff = x[0] - t;
      return Vectord<1>{
        50 * e * diff * (e * diff * diff + 1)
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
