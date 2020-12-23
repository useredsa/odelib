#ifndef INCLUDE_SOLVERS_NEWTOWN_IMPLICIT_SOLVER_HPP_
#define INCLUDE_SOLVERS_NEWTOWN_IMPLICIT_SOLVER_HPP_

#include <algorithm>
#include "initial_value_problem.hpp"
#include "methods/interfaces/plain_implicit_method.hpp"
#include "methods/interfaces/backward_differentiation_formula.hpp"
#include "ode_solution.hpp"
#include "solvers/types.hpp"

namespace odelib {

struct Newton1d {
  static constexpr int kMaxIt = 30;
  static constexpr double kDefTol = 1e-6;

  template <DerivableImplicitEquation Eq>
  inline std::pair<Vectord<1>, bool> solve(const Eq& f, const Vectord<1>& start,
      double tol = kDefTol) {
    Vectord<1> x = start;
    for (int i = 0; i < kMaxIt; ++i) {
      Vectord<1> nxt = x - f(x) / f.dv(x)(0,0);
      if ((x - nxt).norm() < tol) {
        return {nxt, true};
      }
      x = nxt;
    }
    // Newton did not converge
    return {x, false};
  }
};

template <OdeSolution Sol>
bool SuitedForPlainImplicitMethod(const Sol& sol, const SizeArgs& args) {
  if (sol.empty()) {
    std::cerr << "PlainImplicitMethod: you must provide a non-empty solution"
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
SolverResult NewtonAppendNSteps(Sol& sol, const Met& met, const D& f, double h,
    double tol, int n) {
  Newton1d solver;
  double t = sol.t.back();
  auto& x = sol.x;
  int zero = sol.size();
  n += zero;
  sol.resize(n);
  for (int i = zero; i < n; ++i) {
    auto [y, converged] =
        solver.solve(met.equation(f, t, x[i-1], h), x[i-1], tol);
    if (!converged) {
      return SolverResult::kFailedToSolveImplicitEq;
    }
    sol.x[i] = y;
    sol.t[i] = t += h;
    sol.dv[i] = f(sol.t[i], sol.x[i]);
  }
  return SolverResult::kOk;
}

template <IvpDerivative D, PlainImplicitMethod Met, OdeSolution Sol>
SolverResult NewtonExtendPastMaxTime(Sol& sol, const Met& met, const D& f,
    const SizeArgs& args) {
  if (!SuitedForPlainImplicitMethod(sol, args)) {
    return SolverResult::kViolatedPrecondition;
  }
  double h = args.fixedStepSize;
  int iter = std::max(std::ceil((args.maxTime - sol.t.back())/h), 0.0);
  return NewtonAppendNSteps(sol, met, f, h, args.tolerance, iter);
}

template <IvpDerivative D, PlainImplicitMethod Met, OdeSolution Sol,
    CrossFunction Cross>
SolverResult NewtonExtendPastZero(Sol& sol, const Met& met, const D& f,
    const SizeArgs& args, const Cross& cross) {
  if (!SuitedForPlainImplicitMethod(sol, args)) {
    return SolverResult::kViolatedPrecondition;
  }

  Newton1d solver;
  double t = sol.t.back();
  double h = args.fixedStepSize;
  double sgn0 = cross(t, sol.x.back());
  while (t < args.maxTime) {
    auto [y, converged] =
        solver.solve(met.equation(f, t, sol.x.back(), h), sol.x.back(), args.tolerance);
    if (!converged) {
      return SolverResult::kFailedToSolveImplicitEq;
    }
    t += h;
    sol.addPoint(t, y);
    double sgn1 = cross(t, sol.x.back());
    if (sgn0*sgn1 < 0) {
      return SolverResult::kOk;
    }
    sgn0 = sgn1;
  }
  return SolverResult::kExhaustedInterval;
}

template <OdeSolution Sol>
bool SuitedForBdf(const Sol& sol, const SizeArgs& args, int nsteps) {
  if ((int) sol.size() < nsteps + 1) {
    std::cerr << "Bdf: you must provide a solution with " << nsteps
        << " steps!" << std::endl;
    return false;
  }
  if (args.tolerance <= 0) {
    std::cerr << "Bdf: tolerance must be > 0!" << std::endl;
    return false;
  }
  if (args.fixedStepSize <= 0) {
    std::cerr << "Bdf: fixedStepSize must be > 0!" << std::endl;
    return false;
  }
  return true;
}

template <IvpDerivative D, BackwardDifferentiationFormula Met, OdeSolution Sol>
SolverResult NewtonExtendPastMaxTime(Sol& sol, const Met& met, const D& f,
    const SizeArgs& args) {
  constexpr int nsteps = met.kNeededSteps;
  if (!SuitedForBdf(sol, args, nsteps)) {
    return SolverResult::kViolatedPrecondition;
  }
  Newton1d solver;
  double h = args.fixedStepSize;
  double tol = args.tolerance;
  auto& t = sol.t;
  auto& x = sol.x;
  int zero = sol.size();
  int iter = zero + std::max(std::ceil((args.maxTime - t.back())/h), 0.0);
  sol.resize(iter);
  for (int i = zero; i < iter; ++i) {
    auto [y, converged] =
        solver.solve(met.equation(f, t[i-1], &x[i-1]-nsteps, h), x[i-1], tol);
    if (!converged) {
      return SolverResult::kFailedToSolveImplicitEq;
    }
    x[i] = y;
    t[i] = t[i-1] + h;
  }
  return SolverResult::kOk;
}

template <IvpDerivative D, BackwardDifferentiationFormula Met, OdeSolution Sol,
    CrossFunction Cross>
SolverResult NewtonExtendPastZero(Sol& sol, const Met& met, const D& f,
    const SizeArgs& args, const Cross& cross) {
  if (!SuitedForPlainImplicitMethod(sol, args)) {
    return SolverResult::kViolatedPrecondition;
  }

  Newton1d solver;
  double h = args.fixedStepSize;
  auto& t = sol.t;
  auto& x = sol.x;
  auto& dv = sol.dv;
  double sgn0 = cross(t, sol.x.back());
  while (t.back() < args.maxTime) {
    auto [y, converged] =
        solver.solve(met.equation(f, t, x.back(), h), x.back(), args.tolerance);
    if (!converged) {
      return SolverResult::kFailedToSolveImplicitEq;
    }
    t += h;
    sol.addPoint(t, y);
    double sgn1 = cross(t, sol.x.back());
    if (sgn0*sgn1 < 0) {
      return SolverResult::kOk;
    }
    sgn0 = sgn1;
  }
  return SolverResult::kExhaustedInterval;
}


}  // namespace odelib

#endif  // INCLUDE_SOLVERS_NEWTOWN_IMPLICIT_SOLVER_HPP_
