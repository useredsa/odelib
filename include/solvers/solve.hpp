#ifndef INCLUDE_SOLVERS_SOLVE_HPP_
#define INCLUDE_SOLVERS_SOLVE_HPP_

#include "initial_value_problem.hpp"
#include "solutions/standard_ode_solution.hpp"
#include "methods/backwards_euler.hpp"
#include "methods/trapezoidal.hpp"
#include "methods/plain_rk_methods.hpp"
#include "solvers/plain_method_solver.hpp"
#include "solvers/plain_adaptive_method_solver.hpp"
#include "solvers/plain_multistep_method_solver.hpp"
#include "solvers/adaptive_multistep_method_solver.hpp"
#include "solvers/newton_implicit_solver.hpp"
#include "solvers/secant_implicit_solver.hpp"

namespace odelib {

template <typename Met, InitialValueProblem Ivp>
StandardOdeSolution<Ivp::Dv::kDim> SolvePastMaxTime(const Ivp& ivp,
    const Met& met, const SizeArgs& args) {
  StandardOdeSolution sol = StandardOdeSolutionFromIvp(ivp);
  auto result = ExtendPastMaxTime(sol, met, typename Ivp::Dv(), args);
  LogResult(result);
  return sol;
}

template <PlainMultistepMethod Met, InitialValueProblem Ivp>
StandardOdeSolution<Ivp::Dv::kDim> SolvePastMaxTime(const Ivp& ivp,
    const Met& met, const SizeArgs& args) {
  StandardOdeSolution sol = StandardOdeSolutionFromIvp(ivp);
  typename Ivp::Dv f;
  AppendNSteps(sol, RungeKutta<met.kOrder>(), f, args.fixedStepSize,
      met.kNeededSteps);
  auto result = ExtendPastMaxTime(sol, met, f, args);
  LogResult(result);
  return sol;
}

template <AdaptiveMultistepMethod Met, InitialValueProblem Ivp>
StandardOdeSolution<Ivp::Dv::kDim> SolvePastMaxTime(const Ivp& ivp,
    const Met& met, const SizeArgs& args) {
  StandardOdeSolution sol = StandardOdeSolutionFromIvp(ivp);
  RungeKutta<met.kOrder> init;
  auto result = ExtendPastMaxTime(sol, met, init, typename Ivp::Dv(), args);
  LogResult(result);
  return sol;
}

template <PlainImplicitMethod Met, InitialValueProblem Ivp>
StandardOdeSolution<Ivp::Dv::kDim> SolvePastMaxTime(const Ivp& ivp,
    const Met& met, const SizeArgs& args) {
  if constexpr (!SpaceDerivableIvpDerivative<typename Ivp::Dv>) {
    return ForceSecantSolvePastMaxTime(ivp, met, args);
  }
  StandardOdeSolution sol = StandardOdeSolutionFromIvp(ivp);
  typename Ivp::Dv f;
  SolverResult result = NewtonExtendPastMaxTime(sol, met, f, args);
  LogResult(result);
  return sol;
}

template <PlainImplicitMethod Met, InitialValueProblem Ivp>
StandardOdeSolution<Ivp::Dv::kDim> ForceSecantSolvePastMaxTime(const Ivp& ivp,
    const Met& met, const SizeArgs& args) {
  StandardOdeSolution sol = StandardOdeSolutionFromIvp(ivp);
  typename Ivp::Dv f;
  // RungeKutta<met.kOrder> init;
  // AppendNSteps(sol, init, f, args.fixedStepSize, 1);
  double h = args.fixedStepSize;
  double t = sol.t.back();
  auto x = sol.x.back();
  auto c = x.norm() < 1e-10? 1e-18*x.Ones() : (1e-6*x).eval();
  auto [y, converged] = Secant1d().solve(met.equation(f, t, x, h),
      x - c, x, args.tolerance);
  if (!converged) {
    LogResult(SolverResult::kFailedToSolveImplicitEq);
  } else {
    sol.addPoint(t + h, y);
    SolverResult result = SecantExtendPastMaxTime(sol, met, f, args);
    LogResult(result);
  }
  return sol;
}

template <BackwardDifferentiationFormula Met, InitialValueProblem Ivp>
StandardOdeSolution<Ivp::Dv::kDim> SolvePastMaxTime(const Ivp& ivp,
    const Met& met, const SizeArgs& args) {
  if constexpr (!SpaceDerivableIvpDerivative<typename Ivp::Dv>) {
    return ForceSecantSolvePastMaxTime(ivp, met, args);
  }
  StandardOdeSolution sol = StandardOdeSolutionFromIvp(ivp);
  typename Ivp::Dv f;
  SolverResult result = NewtonAppendNSteps(sol, Trapezoidal(), f,
      args.fixedStepSize, args.tolerance, met.kNeededSteps);
  if (result != SolverResult::kOk) {
    LogResult(result);
    return sol;
  }
  // SolverResult result;
  // AppendNSteps(sol, RK4(), f,
  //     args.fixedStepSize, met.kNeededSteps);
  result = NewtonExtendPastMaxTime(sol, met, f, args);
  LogResult(result);
  return sol;
}

template <BackwardDifferentiationFormula Met, InitialValueProblem Ivp>
StandardOdeSolution<Ivp::Dv::kDim> ForceSecantSolvePastMaxTime(const Ivp& ivp,
    const Met& met, const SizeArgs& args) {
  StandardOdeSolution sol = StandardOdeSolutionFromIvp(ivp);
  typename Ivp::Dv f;
  // RungeKutta<met.kOrder> init;
  // AppendNSteps(sol, init, f, args.fixedStepSize, 1);
  double h = args.fixedStepSize;
  double t = sol.t.back();
  auto x = sol.x.back();
  auto c = x.norm() < 1e-10? 1e-18*x.Ones() : (1e-6*x).eval();
  auto [y, converged] = Secant1d().solve(BackwardsEuler().equation(f, t, x, h),
      x - c, x, args.tolerance);
  if (!converged) {
    LogResult(SolverResult::kFailedToSolveImplicitEq);
    return sol;
  }
  sol.addPoint(t + h, y);
  SolverResult result = SecantAppendNSteps(sol, BackwardsEuler(), f, h,
    args.tolerance, met.kNeededSteps - 1);
  if (result != SolverResult::kOk) {
    LogResult(result);
    return sol;
  }
  result = SecantExtendPastMaxTime(sol, met, f, args);
  LogResult(result);
  return sol;
}

}  // namespace odelib

#endif  // INCLUDE_SOLVERS_SOLVE_HPP_

