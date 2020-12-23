#ifndef INCLUDE_PROBLEMS_STOCK_HPP_
#define INCLUDE_PROBLEMS_STOCK_HPP_

#include <cmath>
#include "types.hpp"

namespace odelib {

struct Stock {
  static constexpr double alpha = 0.5;
  static constexpr double gamma = 0.3;
  static constexpr double sigma = 3;
  static constexpr double delta = 0;
  static constexpr double mu = 0.2;
  static constexpr double rho = (1-alpha)*(sigma-1)/(1+alpha*(sigma-1));
  static constexpr double kappa1 = 0.5;
  static constexpr double kappa2 = kappa1;
  static constexpr double kappa = kappa1/kappa2;

  static double phi;

  Stock(const Vectord<3>& x0 = {0.4, 0.8, 0.6}) : t0_(0), x0_(x0) {}

  inline double t0() const { return t0_; }
  inline Vectord<3> x0() const { return x0_; }

  static inline double SolveForW(double L, double tol = 1e-15) {
    double lo = 0, hi = 1e18;
    while (lo + tol < hi) {
      double med = (lo+hi)/2;
      double preval = std::pow(med, alpha*(sigma-1))
          * std::pow(kappa, (1-alpha)*(sigma-1));
      double val = phi*(1-L)*preval + L - med*preval*(phi*L + (1-L)*preval);
      if (val > 0) {
        lo = med;
      } else {
        hi = med;
      }
    }
    return lo;
  }

  struct Dv {
    static constexpr int kDim = 3;

    inline Vectord<3> operator()(double t, const Vectord<3>& x) const {
      constexpr double coef0 = alpha + sigma/(sigma-1);
      constexpr double coef1 = std::pow(kappa, 1 - alpha);
      constexpr double coef2 = rho*(1-kappa1)/kappa1;
      constexpr double coef3 = rho*delta*(1-kappa2)/kappa2;
      constexpr double coef4 = rho*(1-delta)*(1-kappa2)/kappa2;

      double w = SolveForW(x[0]);
      return {
        x[0]*(1-x[0])*(std::pow(w, coef0)*coef1*std::pow(x[1]/x[2], -gamma)-1),
        coef2*w*x[0] + coef3*(1 - x[0]) - mu*x[1],
        coef4*(1 - x[0]) - mu*x[2]
      };
    }
  };

 private:
  double t0_;
  Vectord<3> x0_;
};

}  // namespace odelib

#endif  // INCLUDE_PROBLEMS_STOCK_HPP_
