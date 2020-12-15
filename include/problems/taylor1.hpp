#ifndef INCLUDE_PROBLEMS_TAYLOR_EXAMPLE_HPP_
#define INCLUDE_PROBLEMS_TAYLOR_EXAMPLE_HPP_

#include <cmath>
#include "types.hpp"

namespace odelib {

struct Taylor1 {
  static inline double t0() { return 0; }
  static inline Vectord<1> x0() { return Vectord<1>{0.5}; }

  struct Dv {
    static constexpr int kDim = 1;
    // All derivatives after the 2nd one are equal
    static constexpr int kNDerivableOrder = 1e9;

    inline Vectord<1> operator()(double t, const Vectord<1>& x) const {
      return Vectord<1>{x[0] -t*t +1};
    }

    template<int Order>
    inline Vectord<1> dvn(double t, const Vectord<1>& x) const;
  };
};

template<>
Vectord<1> Taylor1::Dv::dvn<0>(double t, const Vectord<1>& x) const {
  return operator()(t, x);
}

template<>
Vectord<1> Taylor1::Dv::dvn<1>(double t, const Vectord<1>& x) const {
  return Vectord<1>{x[0] - t*t - 2*t + 1};
}

template<int Order>
Vectord<1> Taylor1::Dv::dvn(double t, const Vectord<1>& x) const {
  return Vectord<1>{x[0] - t*t - 2*t - 1};
}

// Alternative Implementation

// template<int Order>
// Vectord<1> Taylor1::Dv::dvn(double t, const Vectord<1>& x) {
//   if constexpr (Order == 0) {
//     return dv(t, x);
//   } else if constexpr (Order == 1) {
//     return Vectord<1>{x[0] - t*t - 2*t + 1};
//   } else {
//     return Vectord<1>{x[0] - t*t - 2*t - 1};
//   }
// }

}  // namespace odelib

#endif  // INCLUDE_PROBLEMS_TAYLOR_EXAMPLE_HPP_
