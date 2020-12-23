#ifndef INCLUDE_METHODS_RICHARDSON_EXTRAPOLATION_HPP_
#define INCLUDE_METHODS_RICHARDSON_EXTRAPOLATION_HPP_

#include <algorithm>
#include <iostream>
#include <utility>
#include <vector>
#include "methods/interfaces/plain_method.hpp"
#include "initial_value_problem.hpp"
#include "types.hpp"

namespace odelib {

/**
 * Richardson Extrapolation for a PlainMethod
 * Increases by one the order of the base method and
 * provides a measure of the error,
 * which makes it an adaptive method.
 * 
 * Computes the next point with step h and h/2 using the base method:
 * (fullstep) x_{n+1} = met.step(x_{n}, h)
 * (halfstep) x_{n+1} = met.step(met.step(x_{n}, h/2), h/2)
 * and extrapolates to a higher order approximation
 * (extrapolation) x_{n+1} = ((2^order)*halfstep - fullstep) / (2^order - 1)
 * giving an approximation of the error as
 * err = (2^order)/(2^order -1) * norm(fullstep - halfstep)
 */
template <PlainMethod Method>
struct RichardsonExtrapolation {
  static constexpr int kOrder = Method::kOrder+1;

  Method met;

  template <IvpDerivative D>
  inline std::pair<Vectord<D::kDim>, double> step(D f, double t,
      const Vectord<D::kDim>& x, double& h, double tolerance) const {
    Vectord<D::kDim> d = f(t, x);
    Vectord<D::kDim> fullstep = met.hinted_step(f, t, x, h, d);
    Vectord<D::kDim> halfstep = met.hinted_step(f, t, x, h/2, d);
    halfstep = met.step(f, t + h/2, halfstep, h/2);

    Vectord<D::kDim> extrapolation = ((1 << Method::kOrder)*halfstep - fullstep)
                               / ((1 << Method::kOrder) - 1);
    constexpr double coef = (1 << Method::kOrder)
                              / static_cast<double>((1 << Method::kOrder) - 1);
    double error = coef*(fullstep-halfstep).norm();

    double q = std::pow(tolerance*h / error, 1.0/Method::kOrder);
    q = std::max(0.1, std::min(q, 4.0));
    h *= q;
    return {extrapolation, error};
  }
};

}  // namespace odelib

#endif  // INCLUDE_METHODS_RICHARDSON_EXTRAPOLATION_HPP_
