#ifndef INCLUDE_METHODS_MULTISTEP_HPP_
#define INCLUDE_METHODS_MULTISTEP_HPP_

#include "initial_value_problem.hpp"
#include "types.hpp"

namespace odelib {

struct AdamsBashforth4 {
  static constexpr int order = 4;
  static constexpr int neededSteps = 3;

  template<IvpDerivative D>
  inline Vectord<D::kDim> step(D f, double t, const Vectord<D::kDim>* x,
      double h, const Vectord<D::kDim>* d) const {
    return x[3] + h/24*(-9*d[0] + 37*d[1] - 59*d[2] + 55*d[3]);
  }
};

}  // namespace odelib

#endif  // INCLUDE_METHODS_MULTISTEP_HPP_
