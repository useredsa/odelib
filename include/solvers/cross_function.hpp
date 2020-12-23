#ifndef INCLUDE_SOLVERS_CROSS_FUNCTION_HPP_
#define INCLUDE_SOLVERS_CROSS_FUNCTION_HPP_

#include <iostream>
#include <utility>
#include "methods/plain_rk_methods.hpp"

namespace odelib {

/**
 * A CrossFunction.
 * A function which defines the a manifold an ODE can cross.
 * The function is positive at one side of the manifold and
 * negative on the other side,
 * being zero inside.
 */
template <typename F>
concept CrossFunction = requires(F f) {
  { F::kDim } -> std::same_as<const int&>;
} && requires(F f, const double t, const Vectord<F::kDim>& x) {
  // Axiom: f is continuous
  { f(t, x) } -> std::same_as<double>;
};

/**
 * Given a point (t0, x0 = x(t0)) where and another time t1 > t0,
 * interval in which the function changes sign,
 * return an approximation via bisection of where (t, x(t)) is 0
 * using a fixed step method.
 */
template <CrossFunction F, IvpDerivative D, PlainMethod Method>
std::pair<double, Vectord<F::kDim>> CrossPoint(const F& cross, const D& f,
    double t0, const Vectord<F::kDim>& x0, double t1, double tol,
    const Method& met) {
  bool sgn = cross(t0, x0) > 0;
  Vectord<F::kDim> dv = f(t0, x0);
  Vectord<F::kDim> x = x0;
  double t = t0, l = t0, r = t1;
  while (l + tol < r) {
    t = (l + r)/2;
    x = met.hinted_step(f, t0, x0, t - t0, dv);
    if ((cross(t, x) > 0) == sgn) {
      l = t;
    } else {
      r = t;
    }
  }
  return {t, x};
}

/**
 * Given a point (t0, x0 = x(t0)) where and another time t1 > t0,
 * interval in which the function changes sign,
 * return an approximation of a given order
 *  via bisection of where (t, x(t)) is 0.
 */
template <int Order, CrossFunction F, IvpDerivative D>
std::pair<double, Vectord<F::kDim>> CrossPoint(const F& cross, const D& f,
    double t0, Vectord<F::kDim> x0, double t1, double tol) {
  return CrossPoint(cross, f, t0, x0, t1, tol, RungeKutta<Order>());
}

/**
 * Given an OdeSolution where the last point crossed the manifold
 * defined by the cross function,
 * returns an approximation of the (time, state) where that happens
 * using a plain method.
 */
template <CrossFunction F, IvpDerivative D, OdeSolution Sol,
    PlainMethod Method>
std::pair<double, Vectord<Sol::kDim>> CrossPoint(const F& cross, const D& f,
    const Sol& sol, double tol, const Method& met) {
  size_t sz = sol.size();
  if (sz < 2) {
    std::cerr << "CrossPoint: The solution contains less than 2 points"
        << std::endl;
  }
  return CrossPoint(cross, f, sol.t[sz-2], sol.x[sz-2], sol.t[sz-1],
      tol, met);
}

/**
 * Given an OdeSolution where the last point crossed the manifold
 * defined by the cross function,
 * returns an approximation of the (time, state) where that happens
 * using a method that assures a given order of approximation.
 */
template <int Order, CrossFunction F, IvpDerivative D, OdeSolution Sol>
std::pair<double, Vectord<Sol::kDim>> CrossPoint(const F& cross, const D& f,
    const Sol& sol, double tol) {
  size_t sz = sol.size();
  if (sz < 2) {
    std::cerr << "CrossPoint: The solution contains less than 2 points"
        << std::endl;
  }
  return CrossPoint<Order>(cross, f, sol.t[sz-2], sol.x[sz-2],
    sol.t[sz-1], tol);
}

template <CrossFunction F, IvpDerivative D, OdeSolution Sol,
    PlainMethod Method>
bool AddCrossPoint(const F& cross, const D& f, const Sol& sol, double tol,
    const Method& met) {
  size_t sz = sol.size();
  if (sz < 2) {
    std::cerr << "CrossPoint: The solution contains less than 2 points"
        << std::endl;
    return false;
  }
  auto [t, x] = CrossPoint(cross, f, sol.t[sz-2], sol.x[sz-2], sol.t[sz-1],
      tol, met);
  sol.addPoint(t, x);
  return true;
}

template <int Order, CrossFunction F, IvpDerivative D, OdeSolution Sol>
bool AddCrossPoint(const F& cross, const D& f, const Sol& sol, double tol) {
  size_t sz = sol.size();
  if (sz < 2) {
    std::cerr << "CrossPoint: The solution contains less than 2 points"
        << std::endl;
    return false;
  }
  auto [t, x] = CrossPoint<Order>(cross, f, sol.t[sz-2], sol.x[sz-2],
      sol.t[sz-1], tol);
  sol.addPoint(t, x);
  return true;
}

}  // namespace odelib

#endif  // INCLUDE_SOLVERS_CROSS_FUNCTION_HPP_

