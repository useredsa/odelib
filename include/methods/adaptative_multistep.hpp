#ifndef INCLUDE_METHODS_ADAPTATIVE_MULTISTEP_HPP_
#define INCLUDE_METHODS_ADAPTATIVE_MULTISTEP_HPP_

#include <iostream>
#include <vector>
#include <utility>
#include <cmath>
#include "Vectord.hpp"

namespace adaptative_multistep {

using Eigen::Vectord;

template<typename Derivative, typename Predictor, typename StartMethod>
struct PredictorCorrector4 {
    Derivative f;
    Predictor predictor;
    StartMethod start_method;
    double tol;
    double h;

    template<int N>
    inline void init(
            std::vector<double>& lt,
            std::vector<Vectord<N>, Eigen::aligned_allocator<Vectord<N>>>& lx,
            std::vector<Vectord<N>, Eigen::aligned_allocator<Vectord<N>>>& ld
            ) {
        ld.push_back(f(lt.back(), lx.back()));
    }

    template<int N>
    inline bool start(
            std::vector<double>& lt,
            std::vector<Vectord<N>, Eigen::aligned_allocator<Vectord<N>>>& lx,
            std::vector<Vectord<N>, Eigen::aligned_allocator<Vectord<N>>>& ld
            ) {
        for(int i=0; i<3; ++i) {
            lx.push_back(lx.back() 
                    + start_method.hinted_step(lt.back(), lx.back(), h, ld.back()));
            lt.push_back(lt.back() + h); // this has to be done here
            ld.push_back(f(lt.back(), lx.back()));
        }
        int ret = step(lt, lx, ld);
        if(ret == 2) {
            lt.erase(lt.end()-3, lt.end());
            lx.erase(lx.end()-3, lx.end());
            ld.erase(ld.end()-3, ld.end());
        }
        return ret;
    }

    template<int N>
    inline int step( // returns 0: ok, 1: ok but restart, 2: failed
            std::vector<double>& lt,
            std::vector<Vectord<N>, Eigen::aligned_allocator<Vectord<N>>>& lx,
            std::vector<Vectord<N>, Eigen::aligned_allocator<Vectord<N>>>& ld
            ) {
        Vectord<N> pred = predictor.step( 
                lx.data() + (lx.size() - 4),
                ld.data() + (ld.size() - 4),
                h
                );
        Vectord<N> dpred = f(lt.back() + h, pred);
        Vectord<N> corr = lx.back() + h/24 * (
                9*dpred + 19*ld[ld.size()-1] - 5*ld[ld.size()-2] + ld[ld.size()-3]
                );
        double diff = (pred - corr).norm();
        double err = diff * 19.0 / 270.0;
        double q = 1.5 * pow(tol * h / diff, 0.25);
        q = std::max(0.1, std::min(q, 4.0));
        if(err < tol * h) {
            lt.push_back(lt.back() + h);
            lx.push_back(corr);
            ld.push_back(f(lt.back(), corr));
            if(err*10 < tol*h) {
                h *= q;
                return 1;
            } 
            return 0;
        } else {
            h *= q;
            return 2;
        }
    }

};

template<typename Method, int N>
std::pair< 
    std::vector<double>,
    std::vector<Vectord<N>, Eigen::aligned_allocator<Vectord<N>>>
    >
adaptative_multistep_ode_solver(
    double t0,
    const Vectord<N>& x0,
    Method met,
    double t_max
) {
    std::vector<double> lt;
    std::vector<Vectord<N>, Eigen::aligned_allocator<Vectord<N>>> lx;
    std::vector<Vectord<N>, Eigen::aligned_allocator<Vectord<N>>> ld;
    lt.push_back(t0);
    lx.push_back(x0);
    met.init(lt, lx, ld);
    while(lt.back() <= t_max) {
        while(met.start(lt, lx, ld));
        while(lt.back() <= t_max && !met.step(lt, lx, ld));
    }
    return {std::move(lt), std::move(lx)};
}

}  // namespace adaptative_multistep

#endif  // INCLUDE_METHODS_ADAPTATIVE_MULTISTEP_HPP_
