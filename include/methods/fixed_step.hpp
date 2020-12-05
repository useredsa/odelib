#ifndef INCLUDE_METHODS_FIXED_STEP_HPP_
#define INCLUDE_METHODS_FIXED_STEP_HPP_

#include <iostream>
#include <vector>
#include "Vectord.hpp"

namespace fixed_step {

using Eigen::Vectord;

template<typename Derivative>
struct Euler {
    Derivative f;

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
        return d*h;
    }
};

template<typename Derivative>
struct ModEuler {
    Derivative f;

    static constexpr int order() {
        return 2;
    }

    template<int N>
    inline Vectord<N> step(double t, const Vectord<N>& x, double h) {
        return hinted_step(t, x, h, f(t, x));
    }

    template<int N>
    inline Vectord<N> hinted_step(double t, const Vectord<N>& x,
                                  double h, const Vectord<N>& d) {
        return (d*h + f(t+h, x+d*h)*h)/2;
    }
};

template<typename Derivative>
struct RK4 {
    Derivative f;

    static constexpr int order() {
        return 4;
    }

    template<int N>
    inline Vectord<N> step(double t, const Vectord<N>& x, double h) {
        return hinted_step(t, x, h, f(t, x));
    }

    template<int N>
    inline Vectord<N> hinted_step(double t, const Vectord<N>& x,
                                  double h, const Vectord<N>& d) {
        Vectord<N> k[4];
        k[0] = d*h;
        k[1] = f(t+h/2, x + k[0]/2)*h;
        k[2] = f(t+h/2, x + k[1]/2)*h;
        k[3] = f(t+h, x + k[2])*h;
        return (k[0] + 2*k[1] + 2*k[2] + k[3])/6;
    }
};

template<typename Derivative, int Order>
struct Taylor {
    Derivative f;

    static_assert(Order <= Derivative::order());

    static constexpr int order() {
        return Order;
    }

    template<int N>
    inline Vectord<N> step(double t, const Vectord<N>& x, double h) {
        Vectord<N> next;
        double coef = 1;
        for(int i=1; i <= Order; ++i) {
            coef *= h / i;
            next += coef * f(t, x, i - 1);
        }
        return next;
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

}  // namespace fixed_step

#endif  // INCLUDE_METHODS_FIXED_STEP_HPP_
