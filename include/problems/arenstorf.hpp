#ifndef INCLUDE_PROBLEMS_ARENSTORF_HPP_
#define INCLUDE_PROBLEMS_ARENSTORF_HPP_

#include <cmath>
#include "Eigen/Dense"

namespace odelib {

struct Arenstorf {
  static constexpr double kMMoon = 0.012277471;
  static constexpr double kMEarth = 1-kMMoon;

  static inline double t0() { return 0; }
  static inline Vectord<4> x0() { return {0.994, 0, 0, -2.001585106}; }

  struct Dv {
    static constexpr int kDim = 4;

    inline Vectord<4> operator()(double t, const Vectord<4>& x) const {
      double D1 = ((x[0]+kMMoon)*(x[0]+kMMoon) + x[1]*x[1]);
      D1 *= std::sqrt(D1);
      double D2 = ((x[0]-kMEarth)*(x[0]-kMEarth) + x[1]*x[1]);
      D2 *= std::sqrt(D2);
      return {
        x[2],
        x[3],
        x[0] + 2*x[3] - kMEarth/D1*(x[0]+kMMoon) - kMMoon/D2*(x[0]-kMEarth),
        x[1] - 2*x[2] - kMEarth/D1*x[1] - kMMoon/D2*x[1]
      };
    }
  };
};

}  // namespace odelib

#endif  // INCLUDE_PROBLEMS_ARENSTORF_HPP_
