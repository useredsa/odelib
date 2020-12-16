#ifndef INCLUDE_METHODS_ADAPTATIVE_MULTISTEP_HPP_
#define INCLUDE_METHODS_ADAPTATIVE_MULTISTEP_HPP_

#include <iostream>
#include <vector>
#include <utility>
#include <cmath>
#include "types.hpp"
#include "methods/interfaces/plain_multistep_method.hpp"
#include "methods/adams_bashforth_4.hpp"

namespace odelib {

template<PlainMultistepMethod Predictor = AdamsBashforth4>
struct PredictorCorrector4 {
  static constexpr int order = std::min(4, Predictor::order);
  static constexpr int neededSteps = std::max(3, Predictor::neededSteps);

  Predictor predictor;

  template<IvpDerivative D>
  inline std::pair<Vectord<D::kDim>, double> step(D f, double t,
      const Vectord<D::kDim>* x, double& h,
      const Vectord<D::kDim>* d, double tol) const {
    Vectord<D::kDim> pred = predictor.step(f, t, x, h, d);
    Vectord<D::kDim> dpred = f(t + h, pred);
    Vectord<D::kDim> corr = x[3] + h/24*(d[1] - 5*d[2] + 19*d[3] + 9*dpred);
    double dist = (pred - corr).norm();
    double error = dist * 19.0 / 270.0;
    double q = 1.5 * pow(tol*h/dist, 0.25);
    q = std::max(0.1, std::min(q, 4.0));
    h *= q;
    return {corr, error};
  }

};

}  // namespace odelib

#endif  // INCLUDE_METHODS_ADAPTATIVE_MULTISTEP_HPP_
