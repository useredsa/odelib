#ifndef INCLUDE_SOLVERS_PLAIN_MULTISTEP_METHOD_SOLVER_HPP_
#define INCLUDE_SOLVERS_PLAIN_MULTISTEP_METHOD_SOLVER_HPP_

#include <algorithm>
#include "initial_value_problem.hpp"
#include "methods/interfaces/plain_multistep_method.hpp"
#include "ode_solution.hpp"
#include "solvers/plain_method_solver.hpp"
#include "solvers/types.hpp"

namespace odelib {

template <OdeSolution Sol>
bool SuitedForMultistepMethod(const Sol& sol, const SizeArgs& args, int steps) {
  if ((int) sol.size() < steps + 1) {
    std::cerr << "PlainMultistepMethod: you must provide a solution with "
        << steps << " steps!" << std::endl;
    return false;
  }
  if (args.fixedStepSize <= 0) {
    std::cerr << "PlainMultistepMethod: fixedStepSize must be > 0!"
        << std::endl;
    return false;
  }
  return true;
}

template <PlainMultistepMethod Met, IvpDerivative D, OdeSolution Sol>
SolverResult ExtendPastMaxTime(Sol& sol, const Met& met, const D& f,
    const SizeArgs& args) {
  constexpr int nsteps = met.kNeededSteps;
  if (!SuitedForMultistepMethod(sol, args, nsteps)) {
    return SolverResult::kViolatedPrecondition;
  }
  double h = args.fixedStepSize;
  auto& t = sol.t;
  auto& x = sol.x;
  auto& dv = sol.dv;
  int zero = sol.size();
  int iter = zero + std::max(std::ceil((args.maxTime - t.back())/h), 0.0);
  sol.resize(iter);
  for (int i = zero; i < iter; ++i) {
    x[i] = met.step(f, t[i-1], &x[i-1]-nsteps, h, &dv[i-1]-nsteps);
    t[i] = t[i-1] + h;
    dv[i] = f(t[i], x[i]);
  }
  return SolverResult::kOk;
}

template <IvpDerivative D, PlainMultistepMethod Met, OdeSolution Sol,
    CrossFunction StopCond>
SolverResult ExtendPastZero(const Met& met, Sol& sol, const D& f,
    const SizeArgs& args, const StopCond& cross) {
  constexpr int nsteps = met.kNeededSteps();
  if (!SuitedForMultistepMethod(sol, args, nsteps)) {
    return SolverResult::kViolatedPrecondition;
  }
  double h = args.fixedStepSize;
  auto& t = sol.t;
  auto& x = sol.x;
  auto& dv = sol.dv;
  double sgn0 = cross(t, sol.x.back());
  while (t.back() < args.maxTime) {
    x.push_back(met.step(f, t, &x.back()-nsteps, h, &dv.back()-nsteps));
    t.push_back(t.back() + h);
    dv.push_back(f(t, x));
    double sgn1 = cross(t, sol.x.back()) > 0;
    if (sgn0*sgn1 < 0) {
      return SolverResult::kOk;
    }
    sgn0 = sgn1;
  }
  return SolverResult::kExhaustedInterval;
}

}  // namespace odelib

#endif  // INCLUDE_SOLVERS_PLAIN_MULTISTEP_METHOD_SOLVER_HPP_

