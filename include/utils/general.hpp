#ifndef INCLUDE_UTILS_GENERAL_HPP_
#define INCLUDE_UTILS_GENERAL_HPP_

#include <iostream>
#include <algorithm>
#include <vector>
#include "Vectord.hpp"

using Eigen::Vectord;

template<int N>
double abs_diff_vectors(
    std::vector<Vectord<N>, Eigen::aligned_allocator<Vectord<N>>> lhs,
    std::vector<Vectord<N>, Eigen::aligned_allocator<Vectord<N>>> rhs) {
    if (lhs.size() != rhs.size()) {
        std::cerr << "abs_diff_vectors: size differ. lhs/rhs: ";
        std::cerr << lhs.size() << "/" << rhs.size() << std::endl;
    }
    int n = std::min(lhs.size(), rhs.size());
    double err = 0;
    for (int i = 0; i < n; ++i) {
        err = std::max(err, (lhs[i]-rhs[i]).norm());
    }
    return err;
}

template<typename AnalyticalSol, int N>
double abs_diff_func(
    std::vector<double> t,
    std::vector<Vectord<N>, Eigen::aligned_allocator<Vectord<N>>> x) {

    AnalyticalSol as;
    if (t.size() != x.size()) {
        std::cerr << "abs_diff_vectors: size differ. t/x: ";
        std::cerr << t.size() << "/" << x.size() << std::endl;
    }
    int n = std::min(t.size(), x.size());
    double err = 0;
    for (int i = 0; i < n; ++i) {
        err = std::max(err, (x[i] - as(t[i])).norm());
    }
    return err;
}

#endif  // INCLUDE_UTILS_GENERAL_HPP_

