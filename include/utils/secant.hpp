#ifndef INCLUDE_UTILS_SECANT_HPP_
#define INCLUDE_UTILS_SECANT_HPP_

#include <iostream>
#include <cmath>

template<typename Function>
struct Secant1d {
    static constexpr int kMaxIt = 30;
    static constexpr double kDefTol = 1e-6;

    Function &f;

    Secant1d(Function &f) : f(f) {}

    inline double solve(double x0, double x1, double tol = kDefTol) {
        double fx0 = f(x0), fx1 = f(x1);
        for (int i = 0; i < kMaxIt; ++i) {
            double x2 = (x0*fx1 - x1*fx0)/(fx1 - fx0);
            if (std::abs(x2 - x1) < tol) {
                return x2;
            }
            x0 = x1, fx0 = fx1;
            x1 = x2, fx1 = f(x2);
        }
        std::cerr << "Secant did not converge." << std::endl;
        return x1;
    }
};

#endif  // INCLUDE_UTILS_SECANT_HPP_
