#ifndef INCLUDE_IMPLICIT_EQUATION_HPP_
#define INCLUDE_IMPLICIT_EQUATION_HPP_

#include "initial_value_problem.hpp"
#include "types.hpp"

namespace odelib {

/**
 * ImplicitEquation
 * An implicit equation is an equation of the form D(x) = 0.
 * We represent such an equation via the function D,
 * which we pass to a root finder algorithm.
 */
template<typename T>
concept ImplicitEquation = requires(T t, const Vectord<1>& x) {
  // { t(x) } -> std::same_as<Vectord<N>>;
  { t({}) };
};

/**
 * DerivableImplicitEquation
 * An implicit equation which is derivable.
 */
template<typename T>
concept DerivableImplicitEquation = ImplicitEquation<T> && requires(T t,
    const Vectord<1>& x) {
  // { t.dv(x) } -> std::same_as<Matrixd<N, N>>;
  { t.dv({}) };
};

/**
 * LfImplicitEquation 
 * Implicit equation D(x) linear in x and f(t, x) for a fixed t.

 * Function of the form D(x) = x - cte0 - cte1*f(cte2, x).
 * It can be used to solve the equation x = cte0 + cte1*f(cte2, x)
 * using numerical methods.
 * When f is SpaceDerivableIvpDerivative then the derivative
 * can also be used.
 */
template<IvpDerivative D>
struct LfImplicitEquation {
  D f;
  double mul, time;
  Vectord<D::kDim> sub;

  LfImplicitEquation(const Vectord<D::kDim>& sub, double mul, double time) :
      sub(sub), mul(mul), time(time) {}

  inline Vectord<D::kDim> operator()(const Vectord<D::kDim>& x) {
    return x - sub - mul*f(time, x);
  }

  template<int Order = 1>
  requires SpaceDerivableIvpDerivative<D> && requires { Order == 1; }
  inline Matrixd<D::kDim, D::kDim> dv(const Vectord<D::kDim>& x) {
    return Matrixd<D::kDim, D::kDim>::Ones() - mul*f.pdvx(time, x);
  }
};

}  // namespace odelib

#endif // INCLUDE_IMPLICIT_EQUATION_HPP_

