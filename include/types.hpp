#ifndef INCLUDE_TYPES_HPP_
#define INCLUDE_TYPES_HPP_

#include <vector>
#include "Eigen/Dense"
#include "Eigen/StdVector"

namespace odelib {

template <int N>
using Vectord = Eigen::Matrix<double, N, 1>;

template <int N>
using VecVectord =
    std::vector<Vectord<N>, Eigen::aligned_allocator<Vectord<N>>>;

template <int N, int M>
using Matrixd = Eigen::Matrix<double, N, M>;

}  // namespace odelib

#endif  // INCLUDE_TYPES_HPP_
