#ifndef INCLUDE_ODE_SOLUTION_HPP_
#define INCLUDE_ODE_SOLUTION_HPP_

#include <algorithm>
#include <concepts>
#include <vector>
#include "types.hpp"
#include "tools/interpolation.hpp"

namespace odelib {

// TODO(edsa): developing a container that does not store
// derivatives and integrating it into the rest of the system
// is still pending.
/**
 * A container to store numerical solutions to ODEs.
 * 
 * The container shall support basic operations,
 * as well as providing a way to determine whether
 * the derivatives at the solution's points need to be calculated.
 */
template <typename T>
concept OdeSolution = requires(const T& ns, int i) {
  { T::kDim } -> std::same_as<const int&>;
  { ns.empty() } -> std::same_as<bool>;
  { ns.size() } -> std::same_as<size_t>;
  { ns.t[i] } -> std::convertible_to<double>;
  { ns.x[i] } -> std::same_as<const Vectord<T::kDim>&>;
  { ns.dv[i] } -> std::same_as<const Vectord<T::kDim>&>;
};

template <OdeSolution Os>
std::optional<Vectord<Os::kDim>> Interpolate(const Os& os, double time,
    int order) {
  int n = os.size();
  int pos = lower_bound(os.t.begin(), os.t.end(), time) - os.t.begin();
  if (pos == n || pos == 0) {
    return {};
  }
  if (order >= 0) {
    const int nPoints = 4;
    // If we don't have 4 points then interpolate with the ones we do have.
    if (n < nPoints) {
      return interpolation::Hermite(os.t, os.x, time);
    }
    // Points we try to use to interpolate if they exist
    // First try with (pos-2 pos-1 pos pos+1)
    // then           (pos-1 pos pos+1 pos+2)...
    const int tries[] = {2, 1, 3, 0};
    for (int id : tries) {
      if (pos >= id && pos+nPoints-id <= n) {
        std::vector<double> t(4);
        std::vector<Vectord<Os::kDim>> x(4);
        for (int i = 0; i < 4; ++i) {
          t[i] = os.t[pos-id+i];
          x[i] = os.x[pos-id+i];
        }
        return interpolation::Hermite(t, x, time);
      }
    }
    // We cannot get here because there are at least nPoints points
  }
  return {};
}

/**
 * Returns the maximum difference between two OdeSolution.
 * 
 * The algorithm interpolates at the time of
 * each of the solution's points.
 */
template <OdeSolution Os1, OdeSolution Os2>
double AbsDiff(Os1 lhs, Os2 rhs) {
  double err = 0;
  for (size_t i = 0; i < lhs.size(); ++i) {
    err = std::max(err, (lhs.x[i] - interpolate(rhs, lhs.t[i])).norm());
  }
  for (size_t i = 0; i < rhs.size(); ++i) {
    err = std::max(err, (rhs.x[i] - interpolate(lhs, rhs.t[i])).norm());
  }
  return err;
}

/**
 * Returns the mean difference between two OdeSolution.
 * 
 * The algorithm interpolates at the time of
 * each of the solution's points.
 */
template <OdeSolution Os1, OdeSolution Os2>
double MeanDiff(Os1 lhs, Os2 rhs) {
  double err = 0;
  for (size_t i = 0; i < lhs.size(); ++i) {
    err += (lhs.x[i] - interpolate(rhs, lhs.t[i])).norm();
  }
  for (size_t i = 0; i < rhs.size(); ++i) {
    err += (rhs.x[i] - interpolate(lhs, rhs.t[i])).norm();
  }
  err /= lhs.size()+rhs.size();
  return err;
}

/**
 * Returns the maximum difference between
 * an OdeSolution and a analytical solution.
 */
template <OdeSolution Os, typename AnalyticalSolution>
double AbsDiff(Os os, AnalyticalSolution as) {
  double err = 0;
  for (size_t i = 0; i < os.size(); ++i) {
    err = std::max(err, (os.x[i]-as(os.t[i])).norm());
  }
  return err;
}

/**
 * Returns the mean difference between
 * an OdeSolution and a analytical solution.
 */
template <OdeSolution Os, typename AnalyticalSolution>
double MeanDiff(Os os, AnalyticalSolution as) {
  double err = 0;
  for (size_t i = 0; i < os.size(); ++i) {
    err += (os.x[i]-as(os.t[i])).norm();
  }
  return err/os.size();
}

}  // namespace odelib

#endif  // INCLUDE_ODE_SOLUTION_HPP_
