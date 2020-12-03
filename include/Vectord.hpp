#ifndef INCLUDE_VECTORD_HPP_
#define INCLUDE_VECTORD_HPP_

#include <Eigen/Dense>
#include <Eigen/StdVector>

namespace Eigen {

template<int N>
using Vectord = Eigen::Matrix<double, N, 1>;

}  // namespace Eigen

#endif // INCLUDE_VECTORD_HPP_

