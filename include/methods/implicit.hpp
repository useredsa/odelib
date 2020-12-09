#ifndef INCLUDE_METHODS_IMPLICIT_HPP_
#define INCLUDE_METHODS_IMPLICIT_HPP_

#include <iostream>
#include <vector>
#include "Vectord.hpp"

namespace implicit {

using Eigen::Vectord;

template<typename Derivative, typename ImplicitCalc>
struct BackwardsEuler {
    Derivative f;
    ImplicitCalc ic; // method that approximates the next derivative

    static constexpr int order() {
        return 1;
    }

    template<int N>
    inline Vectord<N> step(double t, const Vectord<N>& x, double h) {
        return hinted_step(t, x, h, f(t, x));
    }

    template<int N>
    inline Vectord<N> hinted_step(double t, const Vectord<N>& x,
                                  double h, const Vectord<N>& d) {

    }
};

template<typename Method, int N>
std::vector<Vectord<N>, Eigen::aligned_allocator<Vectord<N>>>
fixed_step_ode_solver(
    double t0,
    const Vectord<N>& x0,
    Method met,
    double step = 1e-3,
    int iter = 10000
) {
    std::vector<Vectord<N>, Eigen::aligned_allocator<Vectord<N>>> x(iter+1);
    x[0] = x0;
    for (int i = 0; i < iter; ++i) {
        x[i+1] = x[i] + met.step(t0, x[i], step);
        t0 += step;
    }
    return x;
}

}  // namespace implicit

#endif  // INCLUDE_METHODS_IMPLICIT_HPP_
