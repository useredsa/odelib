#ifndef INCLUDE_PROBLEMS_ARENSTORF_HPP_
#define INCLUDE_PROBLEMS_ARENSTORF_HPP_

#include <cmath>
#include "Eigen/Dense"

namespace arenstorf {

using Eigen::Vector4d;
using std::sqrt;

constexpr int kDim = 4;
constexpr double kMMoon = 0.012277471;
constexpr double kMEarth = 1-kMMoon;

struct derivative {
    inline Vector4d operator()(double t, const Vector4d& x) {
        double D1 = ((x[0]+kMMoon)*(x[0]+kMMoon) + x[1]*x[1]);
        D1 *= std::sqrt(D1);
        double D2 = ((x[0]-kMEarth)*(x[0]-kMEarth) + x[1]*x[1]);
        D2 *= std::sqrt(D2);
        return {
            x[2],
            x[3],
            x[0] + 2*x[3] - kMEarth/D1*(x[0]+kMMoon)
                - kMMoon/D2*(x[0]-kMEarth),
            x[1] - 2*x[2] - kMEarth/D1*x[1] - kMMoon/D2*x[1]
        };
    }
};

}  // namespace arenstorf

#endif  // INCLUDE_PROBLEMS_ARENSTORF_HPP_
