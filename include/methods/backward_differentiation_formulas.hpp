#ifndef INCLUDE_METHODS_BACKWARD_DIFFERENTIATION_FORMULAS_HPP_
#define INCLUDE_METHODS_BACKWARD_DIFFERENTIATION_FORMULAS_HPP_

#include "implicit_equation.hpp"
#include "initial_value_problem.hpp"
#include "types.hpp"

namespace odelib {

template<IvpDerivative D>
struct Bdf2 {
  D f;

  static constexpr int order() {
    return 2;
  }

  static constexpr int neededSteps() {
    return 1;
  }

  inline LfImplicitEquation<D> equation(double t, const Vectord<D::kDim>* x,
                                        double h) {
    return LfImplicitEquation<D>(
        - 1/3.0*x[0] + 4/3.0*x[1], 2/3.0*h, t+h);
  }
};

template<IvpDerivative D>
struct Bdf3 {
  D f;

  static constexpr int order() {
    return 3;
  }

  static constexpr int neededSteps() {
    return 2;
  }

  inline LfImplicitEquation<D> equation(double t, const Vectord<D::kDim>* x,
                                        double h) {
    return LfImplicitEquation<D>(
        2/11.0*x[0] - 9/11.0*x[1] + 18/11.0*x[2], 6/11.0*h, t+h);
  }
};

template<IvpDerivative D>
struct Bdf4 {
  D f;

  static constexpr int order() {
    return 4;
  }

  static constexpr int neededSteps() {
    return 3;
  }

  inline LfImplicitEquation<D> equation(double t, const Vectord<D::kDim>* x,
                                        double h) {
    return LfImplicitEquation<D>(
        - 3/25.0*x[0] + 16/25.0*x[1] - 36/25.0*x[2] + 48/25.0*x[3],
        12/25.0*h,
        t+h);
  }
};

template<IvpDerivative D>
struct Bdf5 {
  D f;

  static constexpr int order() {
    return 5;
  }

  static constexpr int neededSteps() {
    return 4;
  }

  inline LfImplicitEquation<D> equation(double t, const Vectord<D::kDim>* x,
                                        double h) {
    return LfImplicitEquation<D>(12/137.0*x[0] - 75/137.0*x[1]
        + 200/137.0*x[2] - 300/137.0*x[3] + 300/137.0*x[4],
        60/137.0*h,
        t+h);
  }
};

template<IvpDerivative D>
struct Bdf6 {
  D f;

  static constexpr int order() {
    return 6;
  }

  static constexpr int neededSteps() {
    return 5;
  }

  inline LfImplicitEquation<D> equation(double t, const Vectord<D::kDim>* x,
                                        double h) {
    return LfImplicitEquation<D>(- 10/147.0*x[0] + 72/147.0*x[1]
        - 225/147.0*x[2] + 400/147.0*x[3] - 450/147.0*x[4] + 360/147.0*x[5],
        60/147.0*h,
        t+h);
  }
};

}  // namespace odelib

#endif  // INCLUDE_METHODS_BACKWARD_DIFFERENTIATION_FORMULAS_HPP_

