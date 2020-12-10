#include "methods/fixed_step.hpp"
#include "methods/implicit.hpp"
#include "problems/rigid_1d.hpp"
#include "utils/general.hpp"
using namespace implicit;
using namespace fixed_step;
using namespace std;

// ASSIGNMENT 10/12/2020

double t0 = 0;
Vectord<1> x0{-1};
Rigid_1d problem(t0, x0);

double rk4adap_error() {
  return 0;
}

double rkf_error() {
  return 0;
}

// double rk4_error(double step_size, double max_t) {
//   int niters = int(max_t / step_size);
//   vector<double> t(niters+1);
//   for (int i = 0; i < niters+1; ++i)
//     t[i] = step_size*i;
//   auto sol = fixed_step_ode_solver(
//       0,
//       x0,
//       RK4<derivative>(),
//       step_size,
//       niters
//       );
//
//   return abs_diff_func<analytical_sol>(t, sol);
// }

double backwards_euler_error(double step_size, double max_t) {
  int niters = int(max_t / step_size);
  vector<double> t(niters+1);
  for (int i = 0; i < niters+1; ++i)
    t[i] = step_size*i;

  BackwardsEuler<Rigid_1d, 1> method(problem);

  auto sol = method.solve(
      step_size,
      niters
      );

  return abs_diff_func<analytical_sol>(t, sol);
}

double trapezoidal_error(double step_size, double max_t) {
  int niters = int(max_t / step_size);
  vector<double> t(niters+1);
  for (int i = 0; i < niters+1; ++i)
    t[i] = step_size*i;

  Trapezoidal<Rigid_1d, 1> method(problem);

  auto sol = method.solve(
      step_size,
      niters
      );

  return abs_diff_func<analytical_sol>(t, sol);
}

int main() {
  double step_size = 0.000001;
  double t_max = 8;
  cout << trapezoidal_error(step_size, t_max) << endl;
  cout << backwards_euler_error(step_size, t_max) << endl;
  // cout << rk4_error(step_size, t_max) << endl;
}

