#ifndef include/solvers/explicit_solvers_hpp_INCLUDED
#define include/solvers/explicit_solvers_hpp_INCLUDED

#include "initial_value_problem.hpp"
#include "ode_solution.hpp"
#include "interfaces.hpp"

namespace odelib {

struct SolverKwArgs {
  double fixedStepSize;
  double minStepAllowed;
  double maxStepAllowed;
  double tolerance;
};

enum class SolverActionEffect {
  kContinue,
  kSolverStop
};
  
enum class SolverResult {
  kOk,
  kExhaustedInterval,
  kViolatedPrecondition,
  kStopped
};

template<PlainMethod Method>
class PlainMethodSolver {

  PlainMethodSolver(Method met) : met_(met) {}

  template<InitialValueProblem Ivp, AbstractNumericalSolution Ns>
  solve(Ivp ivp, double t1, Ns& ns, SolverKwArgs& kwargs) {
    if (kwargs.fixedStepSize <= 0) {
      std::cerr << "SolverKwArgs: fixedStepSize must be > 0!" << std::endl;
      exit(-1);
    }

    Method met;
    double t = ivp.t0(), step = kwargs.fixedStepSize;
    int iter = std::max(std::ceil((t1 - t)/step), (double) 0);
    Container ns(iter+1, ivp.t0(), step);
    ns.x()[0] = ivp.x0();
    for (int i = 0; i < iter; ++i) {
      ns.x()[i+1] = met.step(t, ns.x()[i], step);
      t += step;
    }
    return ns;
  }

  template<SolverEvent Event>
  solve(Ivp ivp, Event event, SolverKwArgs& kwargs) {
    if (kwargs.fixedStepSize <= 0) {
      std::cerr << "SolverKwArgs: fixedStepSize must be > 0!" << std::endl;
      exit(-1);
    }
    if (kwargs.tolerance < 0) {
      std::cerr << "SolverKwArgs: tolerance must be >= 0!" << std::endl;
      exit(-1);
    }
    double t = ivp.t0();
    double step = kwargs.fixedStepSize;
    double tolerance = kwargs.tolerance;
    double maxTime = kwards.maxTime;
    ns.addPoint(ivp.t0(), ivp.x0());
    while (t < maxTime) {
      ns.addPoint(t + step, met.step(t, ns.x().back(), step);
      t += step;
      if (event.dist(ns) < tolerance) {
        auto action = event.activate(ns);
        if (action == kSolverStop) {
          return kSolverStopped;
        }
      }
    }
    return kSolverFinishedInterval;
  }
};

#endif // include/solvers/explicit_solvers_hpp_INCLUDED

