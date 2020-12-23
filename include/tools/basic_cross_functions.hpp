#ifndef INCLUDE_TOOLS_BASIC_CROSS_FUNCTIONS_HPP_
#define INCLUDE_TOOLS_BASIC_CROSS_FUNCTIONS_HPP_

namespace odelib {

/**
 * A cross function that triggers when an axis is found
 */
template <int Index, int Dim>
struct CrossAxis {
  static constexpr int kDim = Dim;

  double operator()(double t, Vectord<Dim> x) const {
    return x[Index];
  }
};

}  // namespace odelib

#endif  // INCLUDE_TOOLS_BASIC_CROSS_FUNCTIONS_HPP_

