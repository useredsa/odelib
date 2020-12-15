#ifndef INCLUDE_METHODS_TAYLOR_HPP_
#define INCLUDE_METHODS_TAYLOR_HPP_

namespace odelib {

template<NDerivableIvpDerivative D, int Order>
//TODO
// requires ({ 0 < Order && Order <= Derivative::kNDerivableOrder() } -> true);
struct Taylor {
  D f;

  static constexpr int order() {
    return Order;
  }

  inline Vectord<D::kDim> step(double t, const Vectord<D::kDim>& x,
                               double h) const {
      return TaylorsExpansion(t, x, h).compute();
  }

  inline Vectord<D::kDim> hinted_step(double t, const Vectord<D::kDim>& x,
                              double h, const Vectord<D::kDim>& dv) const {
      return h*dv + TaylorsExpansion(t, x, h).template compute<2>();
  }

 private:
  // This class is merely an intent of performing
  // a compile time loop unrolling.
  struct TaylorsExpansion {
      D f;
      Vectord<D::kDim> next = Vectord<D::kDim>::Zero();
      const Vectord<D::kDim>& x;
      double t, h, coef = 1;

      TaylorsExpansion(double t, const Vectord<D::kDim>& x, double h)
          : x(x), t(t), h(h) {}

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

