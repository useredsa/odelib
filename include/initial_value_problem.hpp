#ifndef INCLUDE_INITIAL_VALUE_PROBLEM_HPP_
#define INCLUDE_INITIAL_VALUE_PROBLEM_HPP_

#include <cmath>
#include "Eigen/Dense"

using Eigen::Vectord;

template<typename T, int N>
concept InitialValueProblem = requires(T a, double t, Vectord<N> &x) {
  { a.derivative(t, x) } -> std::same_as<Vectord<N>>;
  { a.get_t0() } -> std::same_as<double>;
  { a.get_x0() } -> std::same_as<Vectord<N>>;
};

template<typename T, int N>
concept ExtendedInitialValueProblem = requires(T a, double t, Vectord<N> &x) {
  { a.dy_derivative(t, x) } -> std::same_as<Vectord<N>>;
  InitialValueProblem<T, N>;
};

#endif  // INCLUDE_INITIAL_VALUE_PROBLEM_HPP_
