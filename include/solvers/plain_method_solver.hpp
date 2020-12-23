#ifndef INCLUDE_SOLVERS_PLAIN_METHOD_SOLVER_HPP_
#define INCLUDE_SOLVERS_PLAIN_METHOD_SOLVER_HPP_

#include <algorithm>
#include "initial_value_problem.hpp"
#include "methods/interfaces/plain_method.hpp"
#include "ode_solution.hpp"
#include "solvers/types.hpp"
#include "solvers/cross_function.hpp"

namespace odelib {

template <OdeSolution Sol>
bool SuitedForPlainMethod(const Sol& sol, const SizeArgs& args) {
  if (sol.empty()) {
    std::cerr << "PlainMethod: you must provide a non-empty solution"
        << std::endl;
    return false;
  }
  if (args.fixedStepSize <= 0) {
    std::cerr << "PlainMethod: fixedStepSize must be > 0!" << std::endl;
    return false;
  }
  return true;
}

template <PlainMethod Met, IvpDerivative D, OdeSolution Sol>
void AppendNSteps(Sol& sol, const Met& met, const D& f, double h, int n) {
    double t = sol.t.back();
    int zero = sol.size();
    n += zero;
    sol.resize(n);
    for (int i = zero; i < n; ++i) {
      sol.x[i] = met.step(f, t, sol.x[i-1], h);
      sol.t[i] = t += h;
      sol.dv[i] = f(sol.t[i], sol.x[i]);
    }
}

template <PlainMethod Met, IvpDerivative D, OdeSolution Sol>
SolverResult ExtendPastMaxTime(Sol& sol, const Met& met, const D& f,
    const SizeArgs& args) {
  if (!SuitedForPlainMethod(sol, args)) {
    return SolverResult::kViolatedPrecondition;
  }
  double h = args.fixedStepSize;
  int iter = std::max(std::ceil((args.maxTime - sol.t.back())/h), 0.0);
  AppendNSteps(sol, met, f, h, iter);
  return SolverResult::kOk;
}

template <PlainMethod Met, IvpDerivative D, OdeSolution Sol,
    CrossFunction StopCond>
SolverResult ExtendPastZero(Sol& sol, const Met& met, const D& f,
    const SizeArgs& args, const StopCond& cross) {
  if (!SuitedForPlainMethod(sol, args)) {
    return SolverResult::kViolatedPrecondition;
  }
  double t = sol.t.back();
  double h = args.fixedStepSize;
  double sgn0 = cross(t, sol.x.back());
  while (t < args.maxTime) {
    sol.addPoint(t + h, met.step(f, t, sol.x.back(), h));
    t += h;
    double sgn1 = cross(t, sol.x.back());
    // If different sign -> we have crossed a zero
    // No matter that cross(t0, x0) == 0,
    // 0 x negative = -0 < 0. See IEEE.
    // Which unfortunately doesn't hold when compiling
    // with --fast-math.
    // Decision is needed.
    //TODO
    if (sgn0*sgn1 < 0 || sgn1 == 0) {
      return SolverResult::kOk;
    }
    sgn0 = sgn1;
  }
  return SolverResult::kExhaustedInterval;
}

}  // namespace odelib

#endif  // INCLUDE_SOLVERS_PLAIN_METHOD_SOLVER_HPP_
