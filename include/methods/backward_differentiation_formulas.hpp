#ifndef INCLUDE_METHODS_BACKWARD_DIFFERENTIATION_FORMULAS_HPP_
#define INCLUDE_METHODS_BACKWARD_DIFFERENTIATION_FORMULAS_HPP_

#include "methods/interfaces/implicit_equation.hpp"
#include "initial_value_problem.hpp"
#include "types.hpp"

namespace odelib {

/**
 * Backward Differentiation Formula of order 2
 * An implicit multistep method.
 */
struct Bdf2 {
  static constexpr int kOrder = 2;
  static constexpr int kNeededSteps = 1;

  template <IvpDerivative D>
  inline LfImplicitEquation<D> equation(const D& f, double t,
      const Vectord<D::kDim>* x, double h) const {
    return LfImplicitEquation<D>(- 1/3.0*x[0] + 4/3.0*x[1], 2/3.0*h, t+h);
  }
};

/**
 * Backward Differentiation Formula of order 3
 * An implicit multistep method.
 */
struct Bdf3 {
  static constexpr int kOrder = 3;
  static constexpr int kNeededSteps = 2;

  template <IvpDerivative D>
  inline LfImplicitEquation<D> equation(const D& f, double t,
      const Vectord<D::kDim>* x, double h) const {
    return LfImplicitEquation<D>(
        2/11.0*x[0] - 9/11.0*x[1] + 18/11.0*x[2], 6/11.0*h, t+h);
  }
};

/**
 * Backward Differentiation Formula of order 4
 * An implicit multistep method.
 */
struct Bdf4 {
  static constexpr int kOrder = 4;
  static constexpr int kNeededSteps = 3;

  template <IvpDerivative D>
  inline LfImplicitEquation<D> equation(const D& f, double t,
      const Vectord<D::kDim>* x, double h) const {
    return LfImplicitEquation<D>(
        - 3/25.0*x[0] + 16/25.0*x[1] - 36/25.0*x[2] + 48/25.0*x[3],
        12/25.0*h,
        t+h);
  }
};

/**
 * Backward Differentiation Formula of order 5
 * An implicit multistep method.
 */
struct Bdf5 {
  static constexpr int kOrder = 5;
  static constexpr int kNeededSteps = 4;

  template <IvpDerivative D>
  inline LfImplicitEquation<D> equation(const D& f, double t,
      const Vectord<D::kDim>* x, double h) const {
    return LfImplicitEquation<D>(
        12/137.0*x[0] - 75/137.0*x[1] + 200/137.0*x[2] - 300/137.0*x[3]
        + 300/137.0*x[4],
        60/137.0*h,
        t+h);
  }
};

/**
 * Backward Differentiation Formula of order 6
 * An implicit multistep method.
 */
struct Bdf6 {
  static constexpr int kOrder = 6;
  static constexpr int kNeededSteps = 5;

  template <IvpDerivative D>
  inline LfImplicitEquation<D> equation(const D& f, double t,
      const Vectord<D::kDim>* x, double h) const {
    return LfImplicitEquation<D>(
        - 10/147.0*x[0] + 72/147.0*x[1] - 225/147.0*x[2] + 400/147.0*x[3]
        - 450/147.0*x[4] + 360/147.0*x[5],
        60/147.0*h,
        t+h);
  }
};

namespace bdfs_unwanted_detail {

template <int Order>
struct SelectedMethods {};

template <>
struct SelectedMethods<2> {
  using Method = Bdf2;
};

template <>
struct SelectedMethods<3> {
  using Method = Bdf3;
};

template <>
struct SelectedMethods<4> {
  using Method = Bdf4;
};

template <>
struct SelectedMethods<5> {
  using Method = Bdf5;
};

template <>
struct SelectedMethods<6> {
  using Method = Bdf6;
};

}  // namespace bdfs_unwanted_detail

/**
 * Alias for the Backward Differentiation Formulas that
 * allows to select one given the order.
 */
template <int Order>
using Bdf = bdfs_unwanted_detail::SelectedMethods<Order>::Method;

}  // namespace odelib

#endif  // INCLUDE_METHODS_BACKWARD_DIFFERENTIATION_FORMULAS_HPP_

