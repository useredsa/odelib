#ifndef INCLUDE_TOOLS_INTERPOLATION_HPP_
#define INCLUDE_TOOLS_INTERPOLATION_HPP_

#include <algorithm>
#include <vector>
#include "types.hpp"

namespace odelib {

namespace interpolation {

/**
 * Interpolates a function from R to R^n for which we know f(x[i]) = y[i]
 * at the point x0.
 */
template <int N>
Vectord<N> Hermite(const std::vector<double>& x,
    const std::vector<Vectord<N>>& y, double x0) {
  int sz = x.size();
  Vectord<N> diff_table[sz][sz];
  std::copy(y.begin(), y.end(), diff_table[0]);

  for (int i = 1; i < sz; ++i) {
    for (int j = i; j < sz; ++j) {
      diff_table[i][j] =
        (diff_table[i-1][j] - diff_table[i-1][j-1]) / (x[j] - x[j-i]);
    }
  }

  Vectord<N> y0{};
  double s = 1;
  for (int i = 0; i < sz; ++i) {
    y0 += s*diff_table[i][sz-1];
    s *= x0 - x[i];
  }
  return y0;
}

}  // namespace interpolation

}  // namespace odelib

#endif  // INCLUDE_TOOLS_INTERPOLATION_HPP_

