#ifndef INCLUDE_METHODS_PLAIN_RK_METHODS_HPP_
#define INCLUDE_METHODS_PLAIN_RK_METHODS_HPP_

#include "methods/euler.hpp"
#include "methods/mod_euler.hpp"
#include "methods/ssprk3.hpp"
#include "methods/rk4.hpp"

namespace odelib {

namespace rungekutta_unwanted_detail {

template <int Order>
struct SelectedMethods {};

template <>
struct SelectedMethods<1> {
  using Method = Euler;
};

template <>
struct SelectedMethods<2> {
  using Method = ModEuler;
};

template <>
struct SelectedMethods<3> {
  using Method = SSPRK3;
};

template <>
struct SelectedMethods<4> {
  using Method = RK4;
};

}  // namespace rungekutta_unwanted_detail

/**
 * A selection of small order plain Runge-Kutta methods.
 * 
 * Useful for selecting a plain method of a given order.
 */
template <int Order>
using RungeKutta = rungekutta_unwanted_detail::SelectedMethods<Order>::Method;

}  // namespace odelib

#endif  // INCLUDE_METHODS_PLAIN_RK_METHODS_HPP_
