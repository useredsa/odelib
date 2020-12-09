#ifndef INCLUDE_UTILS_NEWTON_HPP_
#define INCLUDE_UTILS_NEWTON_HPP_

#include <iostream>
#include <cmath>

template<typename Function, typename Derivative>
struct Newton1d {
    static constexpr int kMaxIt = 30;
    static constexpr double kDefTol = 1e-10;

    Function &f;
    Derivative &der;

    Newton1d(Function &f, Derivative &der)
      : f(f), der(der) {}

    inline double solve(
        double start,
        double tol = kDefTol
        ) {
        double x = start;
        for (int i = 0; i < kMaxIt; ++i) {
            double nxt = x - f(x) / der(x);
            if (std::abs(x - nxt) < tol) {
                return nxt;
            }
            x = nxt;
        }
        std::cerr << "Newton did not converge." << std::endl;
        return x;
    }
};

#endif  // INCLUDE_UTILS_NEWTON_HPP_
