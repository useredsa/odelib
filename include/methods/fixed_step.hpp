#ifndef FIXED_STEP_H
#define FIXED_STEP_H

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/StdVector>

namespace fixed_step {

template<int N>
using Vectord = Eigen::Matrix<double, N, 1>;

template<typename F, int N>
struct Euler {
    F f;

    inline Vectord<N> operator()(double t, const Vectord<N>& x, double step) {
        return f(t, x)*step;
    }
};

template<typename F, int N>
struct ModEuler {
    F f;

    inline Vectord<N> operator()(double t, const Vectord<N>& x, double step) {
        Vectord<N> d = f(t, x)*step;
        return (d + f(t+step, x+d)*step)/2;
    }
};

template<typename F, int N>
struct RK4 {
    F f;

    inline Vectord<N> operator()(double t, const Vectord<N>& x, double step) {
        Vectord<N> k[4];
        k[0] = f(t, x)*step;
        k[1] = f(t+step/2.0, x + k[0]/2.0)*step;
        k[2] = f(t+step/2.0, x + k[1]/2.0)*step;
        k[3] = f(t+step, x + k[2])*step;
        return (k[0] + 2.0*k[1] + 2.0*k[2] + k[3])/6.0;
    }
};

template<typename Method, int N>
std::vector<Vectord<N>, Eigen::aligned_allocator<Vectord<N>>> fixed_step_ode_solver(
    double t0,
    const Vectord<N>& x0,
    Method met,
    double step = 1e-3,
    int iter = 10000
) {
    std::vector<Vectord<N>, Eigen::aligned_allocator<Vectord<N>>> x(iter+1);
    x[0] = x0;
    for (int i = 0; i < iter; ++i) {
        x[i+1] = x[i] + met(t0, x[i], step);
        t0 += step;
    }
    return x;
}

} // namespace fixed_step

#endif // FIXED_STEP_H
