#ifndef INCLUDE_PROBLEMS_TWOBODIES_HPP_
#define INCLUDE_PROBLEMS_TWOBODIES_HPP_

#include <cmath>
#include "types.hpp"

namespace odelib {

class TwoBodies {
 public:
  static constexpr double kMSun = 1988.5;
  static constexpr double kMEarth = 5.9724e-3;
  static constexpr double kGSunEarth = 8.649828e-4;

  TwoBodies(double t0, Vectord<4> x0) : t0_(t0), x0_(x0) {}

  struct Dv {
    static constexpr int kDim = 4;

    double coefficient;

    Dv(double G = kGSunEarth, double M1 = kMSun, double M2 = kMEarth) :
        coefficient(-G * (M1 + M2)) {}

    inline Vectord<4> operator()(double t, const Vectord<4>& x) const {
      double pos_norm3 = pow(x[0]*x[0] + x[1]*x[1], 1.5);
      return {
          x[2],
          x[3],
          coefficient * x[0] / pos_norm3,
          coefficient * x[1] / pos_norm3
      };
    }
  };

  inline double t0() { return t0_; }
  inline Vectord<4> x0() { return x0_; }

 private:
  double t0_;
  Vectord<4> x0_;
};

}  // namespace odelib

#endif  // INCLUDE_PROBLEMS_TWOBODIES_HPP_
