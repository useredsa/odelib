#ifndef INCLUDE_METHODS_NEWTON_HPP_
#define INCLUDE_METHODS_NEWTON_HPP_

#include <iostream>
#include <vector>
#include "Vectord.hpp"

using Eigen::Vectord;

const int MAX_IT = 20;
const int DEFAULT_TOL = 1e-10;

template<typename Function, typename Derivative>
Function f;
Derivative der;

struct Newton {
  inline Vectord<1> (double t, Vectord<1>& start, double tol = DEFAULT_TOL) {
    Vectord<1> x = start;
    for(int i=0; i < MAX_IT; ++i) {
      Vectord<1> nxt = x - f(x) / der(x);
      if((x-nxt).norm() < tol)
        return nxt;
      x = nxt;
    }
    std::cerr << "Newton did not converge." << std::endl;
    return x;
  }
};

#endif  // INCLUDE_METHODS_NEWTON_HPP_
