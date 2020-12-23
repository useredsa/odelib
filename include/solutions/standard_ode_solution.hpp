#ifndef INCLUDE_SOLUTIONS_STANDARD_ODE_SOLUTION_HPP_
#define INCLUDE_SOLUTIONS_STANDARD_ODE_SOLUTION_HPP_

#include <vector>
#include "initial_value_problem.hpp"

namespace odelib {

/**
 * StandardOdeSolution
 * 
 * An OdeSolution that stores time, point and derivative
 * in three independent vectors.
 */
template <int N>
struct StandardOdeSolution {
  static constexpr int kDim = N;
  static constexpr bool kStoresDerivatives = true;

  StandardOdeSolution() {}

  inline size_t size() const { return t.size(); }
  inline bool empty() const { return t.empty(); }
  inline bool containsDerivatives() const { return t.size() == dv.size(); }

  inline void reserve(size_t size) {
    t.reserve(size);
    x.reserve(size);
    dv.reserve(size);
  }

  inline void resize(size_t size) {
    t.resize(size);
    x.resize(size);
    dv.resize(size);
  }

  inline void addPoint(double t, const Vectord<N>& x) {
    this->t.push_back(t);
    this->x.push_back(x);
  }

  inline void addPoint(double t, const Vectord<N>& x, const Vectord<N>& dv) {
    this->t.push_back(t);
    this->x.push_back(x);
    this->dv.push_back(dv);
  }

  std::vector<double> t;
  std::vector<Vectord<N>> x;
  std::vector<Vectord<N>> dv;
};

template <InitialValueProblem Ivp>
StandardOdeSolution<Ivp::Dv::kDim> StandardOdeSolutionFromIvp(Ivp ivp) {
  StandardOdeSolution<Ivp::Dv::kDim> sol;
  sol.addPoint(ivp.t0(), ivp.x0(), typename Ivp::Dv()(ivp.t0(), ivp.x0()));
  return sol;
}

}  // namespace odelib

#endif  // INCLUDE_SOLUTIONS_STANDARD_ODE_SOLUTION_HPP_

