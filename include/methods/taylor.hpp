#ifndef INCLUDE_METHODS_TAYLOR_HPP_
#define INCLUDE_METHODS_TAYLOR_HPP_

namespace odelib {

template<int Order>
//TODO
// requires ({ 0 < Order && Order <= Derivative::kNDerivableOrder() } -> true);
struct Taylor {
  static constexpr int order = Order + 1;

  template<NDerivableIvpDerivative D>
  inline Vectord<D::kDim> step(D f, double t, const Vectord<D::kDim>& x,
      double h) const {
    return TaylorsExpansion(f, t, x, h).compute();
  }

  template<NDerivableIvpDerivative D>
  inline Vectord<D::kDim> hinted_step(D f, double t, const Vectord<D::kDim>& x,
      double h, const Vectord<D::kDim>& dv) const {
    return h*dv + TaylorsExpansion(f, t, x, h).template compute<2>();
  }

 private:
  // This class is merely an intent of performing
  // a compile time loop unrolling.
  template<NDerivableIvpDerivative D>
  struct TaylorsExpansion {
    D f;
    Vectord<D::kDim> next = Vectord<D::kDim>::Zero();
    const Vectord<D::kDim>& x;
    double t, h, coef = 1;

    TaylorsExpansion(D f, double t, const Vectord<D::kDim>& x, double h)
      : f(f), x(x), t(t), h(h) {}

    template<int O = 1>
    inline Vectord<D::kDim>& compute() {
      if constexpr (O == Order) {
        return next;
      } else {
        coef *= h/O;
        next += coef * (f.template f<O+1>(t, x));
        return compute<O+1>();
      }
    }
  };
};

}  // namespace odelib

#endif  // INCLUDE_METHODS_TAYLOR_HPP_

