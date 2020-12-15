#ifndef INCLUDE_INITIAL_VALUE_PROBLEM_HPP_
#define INCLUDE_INITIAL_VALUE_PROBLEM_HPP_

#include <concepts>
#include "types.hpp"

namespace odelib {

// Things that must be always known at compile time (for performance),
// like the dimension of the problem,
// are allowed to be regular constexpr class variables
// Other things like the initial value problem's time and point need not
// be known at compile time for performance,
// and declaring them with getters/setters allows to compute different
// problems depending on user input without the need of recompiling.

template<typename D>
concept IvpDerivative = requires(D f,
    const double t, const Vectord<D::kDim>& x) {
  { D::kDim } -> std::same_as<const int&>;
  // We use a variable here because it's possible
  // (although probably not advisable)
  // to use an object with an evaluation counter.
  { f(t, x) } -> std::same_as<Vectord<D::kDim>>;
};

template<typename D>
concept SpaceDerivableIvpDerivative = IvpDerivative<D> && requires(D f,
    const double t, const Vectord<D::kDim>& x) {
  { f.pdvx(t, x) } -> std::same_as<Matrixd<D::kDim, D::kDim>>;
};

template<typename D>
concept NDerivableIvpDerivative = IvpDerivative<D> && requires(D f,
    const double t, const Vectord<D::kDim>& x) {
  { D::kNDerivableOrder } -> std::same_as<const int&>;
  { f.template dvn<1>(t, x) } -> std::same_as<Vectord<D::kDim>>;
};

template<typename F>
concept NDerivableFunction = requires(F f, const double t) {
  { F::kNDerivableOrder } -> std::same_as<const int&>;
  { f.template dvn<1>(t) };
};

template<typename Ivp>
concept InitialValueProblem = IvpDerivative<typename Ivp::Dv> &&
    requires(Ivp ivp, const double t, const Vectord<Ivp::Dv::kDim>& x) {
  { ivp.t0() } -> std::same_as<double>;
  { ivp.x0() } -> std::same_as<Vectord<Ivp::Dv::kDim>>;
};

}  // namespace odelib

#endif  // INCLUDE_INITIAL_VALUE_PROBLEM_HPP_
