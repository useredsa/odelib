#ifndef INCLUDE_PROBLEMS_TWOBODIES_HPP_
#define INCLUDE_PROBLEMS_TWOBODIES_HPP_

#include <cmath>
#include "Eigen/Dense"

namespace two_bodies {

using Eigen::Vector4d;
using std::sqrt;

constexpr int kDim = 4;
constexpr double kMSun = 1988.5;
constexpr double kMEarth = 5.9724e-3;
constexpr double kGSunEarth = 8.649828e-4;

struct derivative {
    double coefficient;
    
    derivative(
            double G = kGSunEarth,
            double M1 = kMSun,
            double M2 = kMEarth
            ) : coefficient(-G * (M1 + M2)) {}

    inline Vector4d operator()(double t, const Vector4d& x) {
        double pos_norm3 = pow(x[0]*x[0] + x[1]*x[1], 1.5);
        return {
            x[2],
            x[3],
            coefficient * x[0] / pos_norm3,
            coefficient * x[1] / pos_norm3
        };
    }
};

}  // namespace two_bodies

#endif  // INCLUDE_PROBLEMS_TWOBODIES_HPP_
