#ifndef INCLUDE_METHODS_IMPLICIT_HPP_
#define INCLUDE_METHODS_IMPLICIT_HPP_

#include <iostream>
#include <vector>
#include "Vectord.hpp"
#include "utils/newton.hpp"
#include "utils/secant.hpp"
#include "initial_value_problem.hpp"

namespace implicit {

using Eigen::Vectord;

template<typename IVP, int N>
  requires ExtendedInitialValueProblem<IVP, N>
  struct BackwardsEuler {
    
    IVP ivp;
    
    BackwardsEuler(IVP ivp)
      : ivp(ivp) {}

    static constexpr int order() {
      return 1;
    }

    inline Vectord<N> equation(double t, const Vectord<N>& x,
        double h, const Vectord<N>& y) {
      return y - x - h*ivp.derivative(t+h, y);
    }

    inline Vectord<N> hinted_equation(double t, const Vectord<N>& x, double h,
        const Vectord<N>& y, const Vectord<N>& d) {
      return y - x - h*ivp.derivative(t+h, y);
    }

    inline Vectord<N> d_equation(double t, const Vectord<N>& x,
        double h, const Vectord<N>& y) {
      return Vectord<N>::Ones() - h*ivp.dy_derivative(t+h, y);
    }

    std::vector<Vectord<N>, Eigen::aligned_allocator<Vectord<N>>> 
    solve(double step = 1e-3, double iter = 100000){
      return newton_ode_solver(ivp.get_t0(), ivp.get_x0(), *this, step, iter);
    }
  };

template<typename IVP, int N>
  requires ExtendedInitialValueProblem<IVP, N>
  struct Trapezoidal {
    
    IVP ivp;
    
    Trapezoidal(IVP ivp)
      : ivp(ivp) {}

    static constexpr int order() {
      return 1;
    }

    inline Vectord<N> equation(double t, const Vectord<N>& x,
        double h, const Vectord<N>& y) {
      return hinted_equation(t, x, h, y, ivp.derivative(t, x));
    }

    inline Vectord<N> hinted_equation(double t, const Vectord<N>& x, double h,
        const Vectord<N>& y, const Vectord<N>& d) {
      return y - x - h/2*(d + ivp.derivative(t+h, y));
    }

    inline Vectord<N> d_equation(double t, const Vectord<N>& x,
        double h, const Vectord<N>& y) {
      return Vectord<N>::Ones() - h/2*ivp.dy_derivative(t+h, y);
    }

    std::vector<Vectord<N>, Eigen::aligned_allocator<Vectord<N>>> 
    solve(double step = 1e-3, double iter = 100000){
      return newton_ode_solver(ivp.get_t0(), ivp.get_x0(), *this, step, iter);
    }
  };


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
    double t = t0;

    auto f = [&](double w) {
      return met.equation(t, x[i], step, Vectord<1>{w})[0];
    };
    auto der = [&](double w) {
      return met.d_equation(t, x[i], step, Vectord<1>{w})[0];
    };
    Newton1d solver(f, der);

    for (i = 0; i < iter; ++i) {
      x[i+1] = Vectord<1>{solver.solve(x[i][0])};
      t += step;
    }
    return x;
  }

// template<typename Method, int N>
// std::vector<Vectord<N>, Eigen::aligned_allocator<Vectord<N>>>
// secant_ode_solver(
//     double t0,
//     const Vectord<N>& x0,
//     Method met,
//     double step = 1e-3,
//     int iter = 10000
// ) {
//     int i = 0;
//     std::vector<Vectord<N>, Eigen::aligned_allocator<Vectord<N>>> x(iter+1);
//     x[0] = x0;
//
//     auto f = [&](double w) {
//       return met.equation(t0, x[i], step, Vectord<1>{w})[0];
//     };
//     auto der = [&](double w) {
//       return met.d_equation(t0, x[i], step, Vectord<1>{w})[0];
//     };
//     Secant1d solver(f);
//     Newton1d aux(f, der);
//     x[1] = Vectord<1>{aux.solve(x[0][0])};
//
//     for (i = 1; i < iter; ++i) {
//         x[i+1] = Vectord<1>{solver.solve(x[i-1][0], x[i][0])};
//         t0 += step;
//     }
//     return x;
// }

}  // namespace implicit

#endif  // INCLUDE_METHODS_IMPLICIT_HPP_
