#ifndef INCLUDE_PROBLEMS_TAYLOR_EXAMPLE_HPP_
#define INCLUDE_PROBLEMS_TAYLOR_EXAMPLE_HPP_

#include <cmath>
#include "Eigen/Dense"

namespace taylor_example {

using Eigen::Vectord;
using std::sqrt;

constexpr int kDim = 1;

struct derivative {
    static constexpr int order() {
        return 1e6; // all derivatives after the 2nd one are equal
    }
    
    inline Vectord<kDim> operator()(double t, const Vectord<kDim> x, int order = 0) {
        switch(order) {
            case 0:
                return Vectord<1>{x[0] -t*t +1};
            case 1:
                return Vectord<1>{x[0] -t*t -2*t +1};
            default:
                return Vectord<1>{x[0] -t*t -2*t -1};
        }
    }
};

}  // namespace taylor_example

#endif  // INCLUDE_PROBLEMS_TAYLOR_EXAMPLE_HPP_
