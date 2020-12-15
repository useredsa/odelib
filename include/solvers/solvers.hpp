#ifndef INCLUDE_SOLVERS_HPP_
#define INCLUDE_SOLVERS_HPP_

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

template<PlainMethod Method, InitialValueProblem Ivp,
    AbstractNumericalSolution Container = FixedStepNumericalSolution<Ivp::kDim>>
Container solve(Ivp ivp, double t1, SolverKwArgs& kwargs) {
  if (kwargs.fixedStepSize <= 0) {
    std::cerr << "SolverKwArgs: fixedStepSize must be >= 0!" << std::endl;
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

template<AdaptiveSingleStepMethod Method, InitialValueProblem Ivp,
    AbstractNumericalSolution Container = StandardNumericalSolution<Ivp::kDim>>
Container solve(Ivp ivp, double t1, SolverKwArgs& kwargs) {
  if (kwargs.minStepAllowed <= 0) {
    std::cerr << "SolverKwArgs: minStepAllowed must be >= 0!" << std::endl;
    exit(-1);
  }
  if (kwargs.maxStepAllowed <= 0) {
    std::cerr << "SolverKwArgs: maxStepAllowed must be >= 0!" << std::endl;
    exit(-1);
  }
  if (kwargs.tolerance <= 0) {
    std::cerr << "SolverKwArgs: tolerance must be >= 0!" << std::endl;
    exit(-1);
  }

  Method met;
  Container ns;
  double t = ivp.t0();
  auto& x = ns.x();
  double h = kwargs.maxStepAllowed;
  double tol = kwargs.tolerance;

  ns.addPoint(t, ivp.x0());
  while (t < t1) {
    double step = h;
    auto [y, err] = met.step(t, x.back(), h, tol);
    if (err < step*tol) {
      t += step;
      ns.addPoint(t, y);
    }
    if (h < kwargs.minStepAllowed) {
      //TODO find if it's appropiate
      // This check is important
      if (step > kwargs.minStepAllowed) {
        std::cerr << "solve: kita!"
                  << std::endl;
        h = kwargs.minStepAllowed;
      } else {
        std::cerr << "solve: adaptive method went below minStepAllowed!"
                  << std::endl;
        break;
      }
    }
    h = std::min(h, kwargs.maxStepAllowed);
  }
  return ns;
}

template<PlainMultistepMethod Method, PlainMethod Initiator,
    InitialValueProblem Ivp,
    AbstractNumericalSolution Container = StandardNumericalSolution<Ivp::kDim>>
Container solve(Ivp ivp, double t1, SolverKwArgs& kwargs) {
  if (kwargs.fixedStepSize <= 0) {
    std::cerr << "SolverKwArgs: fixedStepSize must be >= 0!" << std::endl;
    exit(-1);
  }

  Method met;
  Initiator init;
  typename Ivp::Dv f;
  double t = ivp.t0(), step = kwargs.fixedStepSize;
  int nsteps = met.neededSteps(), iter = std::max(std::ceil((t1 - t)/step), (double) 0);
  Container ns;
  ns.reserve(iter+1);
  ns.addPoint(t, ivp.x0(), f(t, ivp.x0()));
  for (int i = 0; i < std::min(nsteps, iter); ++i) {
    auto y = init.hinted_step(t, ns.x()[i], step, ns.dv()[i]);
    t += step;
    ns.addPoint(t, y, f(t, y));
  }
  for (int i = nsteps; i < iter; ++i) {
    auto y = met.step(t, &ns.x()[i-nsteps], step, &ns.dv()[i-nsteps]);
    t += step;
    ns.addPoint(t, y, f(t, y));
  }
  return ns;
}

template<AdaptiveMultistepMethod Method, PlainMethod Initiator,
    InitialValueProblem Ivp,
    AbstractNumericalSolution Container = StandardNumericalSolution<Ivp::kDim>>
Container solve(Ivp ivp, double t1, SolverKwArgs& kwargs) {
  if (kwargs.minStepAllowed <= 0) {
    std::cerr << "SolverKwArgs: minStepAllowed must be >= 0!" << std::endl;
    exit(-1);
  }
  if (kwargs.maxStepAllowed <= 0) {
    std::cerr << "SolverKwArgs: maxStepAllowed must be >= 0!" << std::endl;
    exit(-1);
  }
  if (kwargs.tolerance <= 0) {
    std::cerr << "SolverKwArgs: tolerance must be >= 0!" << std::endl;
    exit(-1);
  }


  Method met;
  Initiator init;
  typename Ivp::Dv f;
  double h = std::sqrt(kwargs.minStepAllowed*kwargs.maxStepAllowed);
  double tol = kwargs.tolerance;
  int nsteps = met.neededSteps();
  Container ns;
  ns.addPoint(ivp.t0(), ivp.x0(), f(ivp.t0(), ivp.x0()));
  bool recompute = true;
  auto& t = ns.t();
  auto& x = ns.x();
  auto& dv = ns.dv();

  while (t.back() < t1) {
    if (recompute) {
      for (int i = 0; i < nsteps; ++i) {
        auto y = init.hinted_step(t.back(), x.back(), h, dv.back());
        ns.addPoint(t.back()+h, y, f(t.back()+h, y));
      }
    }
    double step = h;
    auto [y, err] =
        met.step(t.back(), &*(x.end()-(nsteps+1)), h, &*(dv.end()-(nsteps+1)), tol);
    if (err < step*tol) {
      ns.addPoint(t.back()+step, y, f(t.back()+step, y));
      recompute = false;
      // Only if the error is really small adapt the step
      if (err > step*h*0.1) {
        h = step;
      }
    } else if (recompute) {
      t.resize(t.size()-nsteps);
      x.resize(x.size()-nsteps);
      dv.resize(dv.size()-nsteps);
    }
    if (h < kwargs.minStepAllowed) {
      std::cerr << "solve: adaptive method went below minStepAllowed!"
                 << std::endl;
      break;
    }
    h = std::min(h, kwargs.maxStepAllowed);
  }
  return ns;
}

template<typename Method, int N>
std::vector<Vectord<N>, Eigen::aligned_allocator<Vectord<N>>>
newton_ode_solver(
  double t0,
  const Vectord<N>& x0,
  Method met,
  double step = 1e-3,
  int iter = 10000
) {
  int i = 0;
  std::vector<Vectord<N>, Eigen::aligned_allocator<Vectord<N>>> x(iter+1);
  x[0] = x0;

  auto f = [&](double w) {
    return met.equation(t0, x[i], step, Vectord<1>{w})[0];
  };
  auto der = [&](double w) {
    return met.d_equation(t0, x[i], step, Vectord<1>{w})[0];
  };
  Newton1d solver(f, der);

  for (i = 0; i < iter; ++i) {
    x[i+1] = Vectord<1>{solver.solve(x[i][0])};
    t0 += step;
  }
  return x;
}

template<typename Method, int N>
std::vector<Vectord<N>, Eigen::aligned_allocator<Vectord<N>>>
secant_ode_solver(
  double t0,
  const Vectord<N>& x0,
  Method met,
  double step = 1e-3,
  int iter = 10000
) {
  int i = 0;
  std::vector<Vectord<N>, Eigen::aligned_allocator<Vectord<N>>> x(iter+1);
  x[0] = x0;

  auto f = [&](double w) {
    return met.equation(t0, x[i], step, Vectord<1>{w})[0];
  };
  auto der = [&](double w) {
    return met.d_equation(t0, x[i], step, Vectord<1>{w})[0];
  };
  Secant1d solver(f);
  Newton1d aux(f, der);
  x[1] = Vectord<1>{aux.solve(x[0][0])};

  for (i = 1; i < iter; ++i) {
    x[i+1] = Vectord<1>{solver.solve(x[i-1][0], x[i][0])};
    t0 += step;
  }
  return x;
}

}  // namespace odelib

#endif  // INCLUDE_SOLVERS_HPP_

