#ifndef INCLUDE_METHODS_ADAPTIVE_HPP_
#define INCLUDE_METHODS_ADAPTIVE_HPP_

#include <algorithm>
#include <iostream>
#include <utility>
#include "Vectord.hpp"

namespace adaptive {

using Eigen::Vectord;

enum StepAck {
    kStepAccepted,
    kStepRejected,
    kStepWentBelowMin
};

template<int N>
struct adap_return_t {
    Vectord<N> point;
    double time;
    StepAck status;
};

template<typename Derivative, typename Method>
struct RichardsonExtrapolation {
    Derivative f;
    Method met;
    double tolerance;
    double min_step;
    double max_step;
    double current_step;

    RichardsonExtrapolation(double tolerance, double min_step, double max_step)
        : tolerance(tolerance), min_step(min_step), max_step(max_step),
          current_step(max_step) { }

    static constexpr int order() {
        return Method::order()+1;
    }

    template<int N>
    adap_return_t<N> step(double t, const Vectord<N>& x) {
        double& h = current_step;
        Vectord<N> d = f(t, x);
        Vectord<N> fullstep = x + met.hinted_step(t, x, h, d);
        Vectord<N> halfstep = x + met.hinted_step(t, x, h/2, d);
        halfstep += met.step(t+h/2, halfstep, h/2);

        constexpr double coef = (1 << Method::order())
                        / static_cast<double>((1 << Method::order())-1);
        double diff = (fullstep-halfstep).norm();
        double error = coef*diff;
        double q = std::pow(tolerance*h / error, 1.0/Method::order());
        q = std::max(0.1, std::min(q, 4.0));

        if (error < tolerance * h) {
            t += h;
            h *= q;
            Vectord<N> extrapolation =
                ((1 << Method::order())*halfstep - fullstep)
                / ((1 << Method::order()) - 1);
            return {extrapolation, t, kStepAccepted};
        }
        h *= q;
        if (h < min_step) {
            return {halfstep, t, kStepWentBelowMin};
        }
        return {halfstep, t, kStepRejected};
    }
};

template<typename Derivative>
struct RKFelhberg {
    Derivative f;
    double tolerance;
    double min_step;
    double max_step;
    double current_step;

    RKFelhberg(double tolerance, double min_step, double max_step)
        : tolerance(tolerance), min_step(min_step), max_step(max_step),
          current_step(max_step) { }

    constexpr int order() {
        return 4;
    }

    template<int N>
    adap_return_t<N> step(double t, const Vectord<N>& x) {
        double& h = current_step;
        Vectord<N> k[6];

        k[0] = h*f(t, x);
        k[1] = h*f(t +   1/4.0*h, x + 1/4.0*k[0]);
        k[2] = h*f(t +   3/8.0*h, x + (3*k[0] + 9*k[1])/32);
        k[3] = h*f(t + 12/13.0*h, x + (1932*k[0] - 7200*k[1] + 7296*k[2])/2197);
        k[4] = h*f(t +         h, x + 439/216.0*k[0] - 8*k[1] + 3680/513.0*k[2]
                                  - 845/4104.0*k[3]);
        k[5] = h*f(t +   1/2.0*h, x - 8/27.0*k[0] + 2*k[1] - 3544/2565.0*k[2]
                                  + 1859/4104.0*k[3] - 11/40.0*k[4]);

        Vectord<N> rk4 = x + 25/216.0*k[0] + 1408/2565.0*k[2] + 2197/4104.0*k[3]
                       - 1/5.0*k[4];
        Vectord<N> rk5 = x + 16/135.0*k[0] + 6656/12825.0*k[2] + 28561/56430.0*k[3]
                       - 9/50.0*k[4] + 2/55.0*k[5];
        double error = (rk5 - rk4).norm();
        double q = std::pow(tolerance * h / (2*error), 0.25);
        q = std::max(0.1, std::min(q, 4.0));

        if (error < tolerance * h) {
            t += h;
            h *= q;
            return {rk5, t, kStepAccepted};
        }
        h *= q;
        if (h < min_step) {
            return {rk5, t, kStepWentBelowMin};
        }
        return {rk5, t, kStepRejected};
    }
};

template<typename Method, int N>
std::pair<std::vector<double>,
          std::vector<Vectord<N>, Eigen::aligned_allocator<Vectord<N>>>>
adaptive_step_ode_solver(
    double t0,
    const Vectord<N>& x0,
    Method met,
    int iter = 10000
) {
    std::vector<double> t(1, t0);
    std::vector<Vectord<N>, Eigen::aligned_allocator<Vectord<N>>> x(1, x0);
    for (int i = 0; i < iter; ++i) {
        auto r = met.step(t.back(), x.back());
        if (r.status == kStepAccepted) {
            t.push_back(r.time);
            x.push_back(r.point);
        } else if (r.status == kStepWentBelowMin) {
            std::cerr << "We messed up!" << std::endl;
            std::cerr << "Iteration: " << i << std::endl;
            break;
        }
    }
    return {std::move(t), std::move(x)};
}

}  // namespace adaptive

#endif  // INCLUDE_METHODS_ADAPTIVE_HPP_
