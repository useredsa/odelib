#ifndef INCLUDE_TOOLS_INTERPOLATION_HPP_
#define INCLUDE_TOOLS_INTERPOLATION_HPP_

#include <algorithm>
#include <vector>
#include "Vecotrd"

namespace interpolation {

using std::vector;
using Eigen::Vecotrd;

template<int N>
hermite(const vector<double>& x, const vector<Vectord<N>>& y,
        const Vectord<N>& x0) {
    int sz = x.size();
    Vectord<N> diff_table[sz][sz];
    std::copy(y.begin(), y.end(), diff_table[0]);

    for (int i = 1; i < sz; ++i) {
        for (int j = i; j < sz; ++j) {
            diff_table[i][j] =
                (diff_table[i-1][j] - diff_table[i-1][j-1]) / (x[j] - x[j-i]);
        }
    }

    Vecotrd<N> y0{};
    double s = 1;
    for (int i = 0; i < sz; ++i) {
        y0 += s*diff_table[i][sz-1];
        s *= x0 - x[i];
    }
    return y0;
}

}  // namespace interpolation

#endif  // INCLUDE_TOOLS_INTERPOLATION_HPP_

