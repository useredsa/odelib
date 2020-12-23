#ifndef INCLUDE_SOLVERS_PLAIN_ADAPTIVE_METHOD_SOLVER_HPP_
#define INCLUDE_SOLVERS_PLAIN_ADAPTIVE_METHOD_SOLVER_HPP_

#include <algorithm>
#include "initial_value_problem.hpp"
#include "methods/interfaces/plain_adaptive_method.hpp"
#include "ode_solution.hpp"
#include "solvers/types.hpp"

namespace odelib {

template <OdeSolution Sol>
bool SuitedForAdaptiveMethod(const Sol& sol, const SizeArgs& args) {
  if (sol.empty()) {
    std::cerr << "PlainAdaptiveMethod: you must provide a non-empty solution"
        << std::endl;
    return false;
  }
  if (args.tolerance <= 0) {
    std::cerr << "PlainAdaptiveMethod: tolerance must be > 0!" << std::endl;
    return false;
  }
  if (args.minStepAllowed <= 0) {
    std::cerr << "PlainAdaptiveMethod: minStepAllowed must be > 0!"
        << std::endl;
    return false;
  }
  if (args.maxStepAllowed <= 0) {
    std::cerr << "PlainAdaptiveMethod: maxStepAllowed must be > 0!"
        << std::endl;
    return false;
  }
  return true;
}

template <IvpDerivative D, PlainAdaptiveMethod Met, OdeSolution Sol>
SolverResult ExtendPastMaxTime(Sol& sol, const Met& met, const D& f,
    const SizeArgs& args) {
  if (!SuitedForAdaptiveMethod(sol, args)) {
    return SolverResult::kViolatedPrecondition;
  }

  // Theoritically, we could start with the minimum step possible
  // because it would automatically increase and we would not fail
  // a few first steps. But according to the solutions obtained
  // for arenstorf's problem, it seems best to start with maximum step.
  // Further investigation is required.
  double h = args.maxStepAllowed;
  double tol = args.tolerance;
  double t = sol.t.back();
  const auto& x = sol.x;
  while (t < args.maxTime) {
    // h is passed by referenced to the method
    // this is the value used in the current step.
    double step = h;
    auto [y, err] = met.step(f, t, x.back(), h, tol);
    if (err < step*tol) {
      t += step;
      sol.addPoint(t, y);
    }
    if (h < args.minStepAllowed) {
      //TODO find if this check is important (and comment here after)
      if (step > args.minStepAllowed) {
        h = args.minStepAllowed;
      } else {
        return SolverResult::kStepWentBelowMin;
      }
    }
    h = std::min(h, args.maxStepAllowed);
  }
  return SolverResult::kOk;
}

template <IvpDerivative D, PlainAdaptiveMethod Met, OdeSolution Sol,
    CrossFunction Cross>
SolverResult ExtendPastMaxTime(Sol& sol, const Met& met, const D& f,
    const SizeArgs& args) {
  if (!SuitedForAdaptiveMethod(sol, args)) {
    return SolverResult::kViolatedPrecondition;
  }

  double h = args.minStepAllowed;
  double tol = args.tolerance;
  double t = sol.t.back();
  const auto& x = sol.x;
  double sgn0 = cross(t, sol.x.back());
  while (t < args.maxTime) {
    double step = h;
    auto [y, err] = met.step(f, t, x.back(), h, tol);
    if (err < step*tol) {
      t += step;
      sol.addPoint(t, y);
      double sgn1 = cross(t, y);
      if (sgn0*sgn1 < 0) {
        return SolverResult::kOk;
      }
      sgn0 = sgn1;
    }
    if (h < args.minStepAllowed) {
      //TODO find if this check is important (and comment here after)
      if (step > args.minStepAllowed) {
        h = args.minStepAllowed;
      } else {
        return SolverResult::kStepWentBelowMin;
      }
    }
    h = std::min(h, args.maxStepAllowed);
  }
  return SolverResult::kExhaustedInterval;
}

}  // namespace odelib

#endif // INCLUDE_SOLVERS_PLAIN_ADAPTIVE_METHOD_SOLVER_HPP_

