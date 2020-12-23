#ifndef INCLUDE_SOLVERS_ADAPTIVE_MULTISTEP_METHOD_SOLVER_HPP_
#define INCLUDE_SOLVERS_ADAPTIVE_MULTISTEP_METHOD_SOLVER_HPP_

#include <algorithm>
#include "initial_value_problem.hpp"
#include "methods/interfaces/adaptive_multistep_method.hpp"
#include "methods/interfaces/plain_method.hpp"
#include "ode_solution.hpp"
#include "solvers/types.hpp"
#include "solvers/plain_method_solver.hpp"

namespace odelib {

template <OdeSolution Sol>
bool SuitedForAdaptiveMultistepMethod(const Sol& sol, const SizeArgs& args,
    int steps) {
  if (static_cast<int>(sol.size()) < steps + 1) {
    std::cerr << "AdaptiveMultistepMethod: you must provide a solution with "
        << steps << " steps!" << std::endl;
    return false;
  }
  if (args.tolerance <= 0) {
    std::cerr << "AdaptiveMultistepMethod: tolerance must be > 0!" << std::endl;
    return false;
  }
  if (args.minStepAllowed <= 0) {
    std::cerr << "AdaptiveMultistepMethod: minStepAllowed must be > 0!"
        << std::endl;
    return false;
  }
  if (args.maxStepAllowed <= 0) {
    std::cerr << "AdaptiveMultistepMethod: maxStepAllowed must be > 0!"
        << std::endl;
    return false;
  }
  return true;
}

template <IvpDerivative D, AdaptiveMultistepMethod Met, PlainMethod Init,
    OdeSolution Sol>
SolverResult ExtendPastMaxTime(Sol& sol, const Met& met, const Init& init,
    const D& f, const SizeArgs& args, bool recompute = true) {
  constexpr int nsteps = met.kNeededSteps;
  if (!SuitedForAdaptiveMultistepMethod(sol, args, recompute? 0 : nsteps)) {
    return SolverResult::kViolatedPrecondition;
  }
  double tol = args.tolerance;
  auto& t = sol.t;
  auto& x = sol.x;
  auto& dv = sol.dv;
  double h = recompute? std::sqrt(args.minStepAllowed*args.maxStepAllowed) :
      t[t.size()-1] - t[t.size()-2];

  while (t.back() < args.maxTime) {
    if (recompute) {
      AppendNSteps(sol, init, f, h, nsteps);
    }
    double step = h;
    auto [y, err] = met.step(f, t.back(), &*(x.end()-(nsteps+1)), h,
        &*(dv.end()-(nsteps+1)), tol);
    if (err < step*tol) {
      x.emplace_back(y);
      t.push_back(t.back()+step);
      dv.push_back(f(t.back(), y));
      recompute = false;
      // If the error is not really small, keep the last step size
      if (err > step*tol*0.1) {
        h = step;
      } else {
        recompute = true;
      }
    } else {
      if (recompute) {
        // Remove the extra steps on failure
        t.resize(t.size()-nsteps);
        x.resize(x.size()-nsteps);
        dv.resize(dv.size()-nsteps);
      }
      recompute = true;
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

template <IvpDerivative D, AdaptiveMultistepMethod Met, PlainMethod Init,
    OdeSolution Sol, CrossFunction StopCond>
SolverResult ExtendPastZero(Sol& sol, const Met& met, const Init& init,
    const D& f, const SizeArgs& args, const StopCond& cross,
    bool recompute = true) {
  constexpr int nsteps = met.kNeededSteps;
  if (!SuitedForAdaptiveMultistepMethod(sol, args, recompute? 0 : nsteps)) {
    return SolverResult::kViolatedPrecondition;
  }

  double tol = args.tolerance;
  auto& t = sol.t();
  auto& x = sol.x();
  auto& dv = sol.dv();
  double h = recompute? std::sqrt(args.minStepAllowed*args.maxStepAllowed) :
      t[t.size()-1] - t[t.size()-2];
  double sgn0 = cross(t.back(), x.back());

  while (t.back() < args.maxTime) {
    if (recompute) {
      stepN(sol, f, init, h, nsteps);
    }
    double step = h;
    auto [y, err] = met.step(t.back(), &*(x.end()-(nsteps+1)), h,
        &*(dv.end()-(nsteps+1)), tol);

    if (err < step*tol) {
      x.emplace_back(y);
      t.push_back(t.back()+step);
      dv.push_back(f(t.back(), y));
      if (recompute) {
        // Check if the cross happened during the extra steps
        for (int i = nsteps+1; i > 0; --i) {
          double sgn1 = cross(t[t.size()-i], x[x.size()-i]);
          if (sgn0*sgn1 < 0) {
            t.resize(t.size()-i+1);
            x.resize(x.size()-i+1);
            dv.resize(dv.size()-i+1);
            return SolverResult::kOk;
          }
          sgn0 = sgn1;
        }
        recompute = false;
      }
      double sgn1 = cross(t.back(), x.back());
      if (sgn0*sgn1 < 0) {
        return SolverResult::kOk;
      }
      sgn0 = sgn1;
      // If the error is not really small, keep the last step size
      if (err > step*tol*0.1) {
        h = step;
      } else {
        recompute = true;
      }
    } else {
      if (recompute) {
        // Remove the extra steps on failure
        t.resize(t.size()-nsteps);
        x.resize(x.size()-nsteps);
        dv.resize(dv.size()-nsteps);
      }
      recompute = true;
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

}  // namespace odelib

#endif  // INCLUDE_SOLVERS_ADAPTIVE_MULTISTEP_METHOD_SOLVER_HPP_
