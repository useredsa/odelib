#ifndef INCLUDE_METHODS_RICHARDSON_EXTRAPOLATION_HPP_
#define INCLUDE_METHODS_RICHARDSON_EXTRAPOLATION_HPP_

#include <algorithm>
#include <iostream>
#include <utility>
#include <vector>
#include "methods/interfaces/adaptive_single_step_method.hpp"
#include "methods/interfaces/plain_method.hpp"
#include "initial_value_problem.hpp"
#include "types.hpp"

namespace odelib {

template<PlainMethod Method, IvpDerivative D>
struct RichardsonExtrapolation {
  Method met;
  D f;

  static constexpr int order() {
    return Method::order()+1;
  }

  inline std::pair<Vectord<D::kDim>, double> step(double t,
      const Vectord<D::kDim>& x, double& h, double tolerance) const {
    Vectord<D::kDim> d = f(t, x);
    Vectord<D::kDim> fullstep = met.hinted_step(t, x, h, d);
    Vectord<D::kDim> halfstep = met.hinted_step(t, x, h/2, d);
    halfstep = met.step(t+h/2, halfstep, h/2);
    Vectord<D::kDim> extrapolation = ((1 << Method::order())*halfstep - fullstep)
                               / ((1 << Method::order()) - 1);

    constexpr double coef = (1 << Method::order())
                              / static_cast<double>((1 << Method::order())-1);
    double error = coef*(fullstep-halfstep).norm();
    double q = std::pow(tolerance*h / error, 1.0/met.order());
    q = std::max(0.1, std::min(q, 4.0));
    h *= q;
    return {extrapolation, error};
  }
};

}  // namespace odelib

#endif  // INCLUDE_METHODS_RICHARDSON_EXTRAPOLATION_HPP_
