#ifndef INCLUDE_SOLVERS_SECANT_IMPLICIT_SOLVER_HPP_
#define INCLUDE_SOLVERS_SECANT_IMPLICIT_SOLVER_HPP_

#include <algorithm>
#include <utility>
#include "initial_value_problem.hpp"
#include "methods/interfaces/plain_implicit_method.hpp"
#include "ode_solution.hpp"
#include "solvers/types.hpp"

namespace odelib {

struct Secant1d {
  static constexpr int kMaxIt = 30;
  static constexpr double kDefTol = 1e-6;

  template <ImplicitEquation Eq>
  inline std::pair<Vectord<1>, bool> solve(Eq f, Vectord<1> x0, Vectord<1> x1,
      double tol = kDefTol) {
    Vectord<1> fx0 = f(x0), fx1 = f(x1);
    for (int i = 0; i < kMaxIt; ++i) {
      Vectord<1> x2 = (x0*fx1 - x1*fx0)/(fx1 - fx0)[0];
      if ((x2 - x1).norm() < tol) {
        return {x2, true};
      }
      x0 = x1, fx0 = fx1;
      x1 = x2, fx1 = f(x2);
    }

    // Secant did not converge
    return {x1, false};
  }
};

// We need two points to use the secant solver
template <OdeSolution Sol>
bool SuitedForSecantPlainImplicitMethod(const Sol& sol, const SizeArgs& args) {
  if ((int) sol.size() < 2) {
    std::cerr << "Secant: you must provide a solution with at least 2 points"
        << std::endl;
    return false;
  }
  if (args.fixedStepSize <= 0) {
    std::cerr << "PlainImplicitMethod: fixedStepSize must be > 0!" << std::endl;
    return false;
  }
  if (args.tolerance <= 0) {
    std::cerr << "PlainImplicitMethod: tolerance must be > 0!" << std::endl;
    return false;
  }
  return true;
}

template <PlainImplicitMethod Met, IvpDerivative D, OdeSolution Sol>
SolverResult SecantAppendNSteps(Sol& sol, const Met& met, const D& f, double h,
    double tol, int n) {
  Secant1d solver;
  double t = sol.t.back();
  auto& x = sol.x;
  int zero = sol.size();
  n += zero;
  sol.resize(n);
  for (int i = zero; i < n; ++i) {
    auto [y, converged] =
        solver.solve(met.equation(f, t, x[i-1], h), x[i-2], x[i-1], tol);
    if (!converged) {
      return SolverResult::kFailedToSolveImplicitEq;
    }
    x[i] = y;
    sol.t[i] = t += h;
    sol.dv[i] = f(sol.t[i], sol.x[i]);
  }
  return SolverResult::kOk;
}

template <IvpDerivative D, PlainImplicitMethod Met, OdeSolution Sol>
SolverResult SecantExtendPastMaxTime(Sol& sol, const Met& met, const D& f,
    const SizeArgs& args) {
  if (!SuitedForSecantPlainImplicitMethod(sol, args)) {
    return SolverResult::kViolatedPrecondition;
  }
  double h = args.fixedStepSize;
  int iter = std::max(std::ceil((args.maxTime - sol.t.back())/h), 0.0);
  return SecantAppendNSteps(sol, met, f, h, args.tolerance, iter);
}

template <IvpDerivative D, PlainImplicitMethod Met, OdeSolution Sol,
    CrossFunction Cross>
SolverResult SecantExtendPastZero(Sol& sol, const Met& met, const D& f,
    const SizeArgs& args, const Cross& cross) {
  if (!SuitedForSecantPlainImplicitMethod(sol, args)) {
    return SolverResult::kViolatedPrecondition;
  }

  Secant1d solver;
  double h = args.fixedStepSize;
  double tol = args.tolerance;
  double maxTime = args.maxTime;
  double t = sol.t.back();
  auto& x = sol.x;
  double sgn0 = cross(t, x.back());
  while (t < maxTime) {
    auto [y, converged] = solver.solve(met.equation(f, t, x.back(), h),
      *(x.end()-2), x.back(), tol);
    if (!converged) {
      return SolverResult::kFailedToSolveImplicitEq;
    }
    t += h;
    sol.addPoint(t, y);
    double sgn1 = cross(t, x.back());
    if (sgn0*sgn1 < 0) {
      return SolverResult::kOk;
    }
    sgn0 = sgn1;
  }
  return SolverResult::kExhaustedInterval;
}

}  // namespace odelib

#endif  // INCLUDE_SOLVERS_SECANT_IMPLICIT_SOLVER_HPP_

