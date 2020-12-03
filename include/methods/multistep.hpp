#ifndef INCLUDE_METHODS_MULTISTEP_HPP_
#define INCLUDE_METHODS_MULTISTEP_HPP_

#include "Vectord.hpp"

namespace multistep {

using Eigen::Vectord;

template<typename Derivative>
struct AdamsBashforth4 {
    Derivative f;

    constexpr int order() {
        return 4;
    }

    constexpr int needed_steps() {
        return 3;
    }

    template<int N>
    inline Vectord<N> step(const Vectord<N>* x, const Vectord<N>* d, double h) {
        return x[3] + h/24*(-9*d[0] + 37*d[1] - 59*d[2] + 55*d[3]);
    }
};

}  // namespace multistep

#endif // INCLUDE_METHODS_MULTISTEP_HPP_

