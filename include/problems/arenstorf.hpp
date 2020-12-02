#ifndef ARENSTORF_H
#define ARENSTORF_H

#include <cmath>
#include <Eigen/Dense>

namespace arenstorf {

using Eigen::Vector4d;
using std::sqrt;

const int K_DIM = 4;
const double K_MMOON = 0.012277471;
const double K_MEARTH = 1-K_MMOON;

struct derivative {
    inline Vector4d operator()(double t, const Vector4d& x) {
        double D1 = ((x[0]+K_MMOON)*(x[0]+K_MMOON) + x[1]*x[1]);
        D1 *= std::sqrt(D1);
        double D2 = ((x[0]-K_MEARTH)*(x[0]-K_MEARTH) + x[1]*x[1]);
        D2 *= std::sqrt(D2);
        return {
            x[2],
            x[3],
            x[0] + 2*x[3] - K_MEARTH/D1*(x[0]+K_MMOON) - K_MMOON/D2*(x[0]-K_MEARTH),
            x[1] - 2*x[2] - K_MEARTH/D1*x[1] - K_MMOON/D2*x[1]
        };
    }
};

} // namespace arenstorf

#endif
