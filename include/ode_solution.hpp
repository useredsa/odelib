#ifndef INCLUDE_ODE_SOLUTION_HPP_
#define INCLUDE_ODE_SOLUTION_HPP_

#include <concepts>
#include <vector>
#include "types.hpp"

namespace odelib {

template<typename T>
concept AbstractNumericalSolution = requires(const T& ns, int i) {
  { T::kDim } -> std::same_as<const int&>;
  { ns.size() } -> std::same_as<size_t>;
  { ns.t()[i] } -> std::convertible_to<double>;
  { ns.x()[i] } -> std::same_as<const Vectord<T::kDim>&>;
  { ns.dv()[i] } -> std::same_as<const Vectord<T::kDim>&>;
};

template<typename T>
concept SpecifiesTolerance = requires(T t) {
  { t.tolerance() } -> std::same_as<const double&>;
};

template<typename T>
concept SupportsFixedStep = requires(T ns) {
  { ns.fixedStepSize() } -> std::same_as<double>;
};

template<typename T>
concept SupportsAdaptableStep = SpecifiesTolerance<T> && requires(T ns) {
  { ns.minStepAllowed() } -> std::same_as<double>;
  { ns.maxStepAllowed() } -> std::same_as<double>;
};

template<typename T>
concept SupportsImplicitSteps = SpecifiesTolerance<T>;

template<int N>
struct StandardNumericalSolution {
 public:
  static constexpr int kDim = N;

  StandardNumericalSolution() {}

  inline size_t size() const { return t_.size(); }
  inline const std::vector<double>& t() const { return t_; }
  inline const std::vector<Vectord<N>>& x() const { return x_; }
  inline const std::vector<Vectord<N>>& dv() const { return dv_; }
  inline std::vector<double>& t() { return t_; }
  inline std::vector<Vectord<N>>& x() { return x_; }
  inline std::vector<Vectord<N>>& dv() { return dv_; }

  inline void reserve(size_t size) {
    t_.reserve(size);
    x_.reserve(size);
    dv_.reserve(size);
  }

  inline void addPoint(double t, const Vectord<N>& x) {
    t_.push_back(t);
    x_.push_back(x);
  }
  inline void addPoint(double t, const Vectord<N>& x, const Vectord<N>& dv) {
    t_.push_back(t);
    x_.push_back(x);
    dv_.push_back(dv);
  }

 private:
  std::vector<double> t_;
  std::vector<Vectord<N>> x_;
  std::vector<Vectord<N>> dv_;
};

struct LinearRange {
 public:
  LinearRange(double t0, double h) : t0_(t0), h_(h) {}

  inline double operator[](int i) const { return t0_ + i*h_; }

 private:
  double t0_, h_;
};

template<int N>
struct FixedStepNumericalSolution {
 public:
  static constexpr int kDim = N;

  FixedStepNumericalSolution(double t0, double h) : t_(t0, h) {}

  FixedStepNumericalSolution(int size, double t0, double h) : t_(t0, h),
      x_(size), dv_(size) {}

  inline size_t size() const { return x_.size(); }
  inline const LinearRange& t() const { return t_; }
  inline const std::vector<Vectord<N>>& x() const { return x_; }
  inline const std::vector<Vectord<N>>& dv() const { return dv_; }
  inline std::vector<Vectord<N>>& x() { return x_; }
  inline std::vector<Vectord<N>>& dv() { return dv_; }

  inline void reserve(size_t size) {
    x_.reserve(size);
    dv_.reserve(size);
  }

  inline void addPoint(double t, const Vectord<N>& x, const Vectord<N>& dv) {
    x_.push_back(x);
    dv_.push_back(dv);
  }

 private:
  LinearRange t_;
  std::vector<Vectord<N>> x_;
  std::vector<Vectord<N>> dv_;
};

}  // namespace odelib

#endif  // INCLUDE_ODE_SOLUTION_HPP_

