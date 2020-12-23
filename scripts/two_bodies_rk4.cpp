#include <iostream>
#include <iomanip>
// Select a problem you want to solve
#include "problems/two_bodies.hpp"
// Include the method you will use in the problem
#include "methods/rk4.hpp"
#include "solvers/plain_method_solver.hpp"
// Include a container for the solution
#include "solutions/standard_ode_solution.hpp"
// Include additional tools
#include "tools/basic_cross_functions.hpp"
#include "tools/tsv_output.hpp"
using namespace std;
using namespace odelib;
using Dv = TwoBodies::Dv;

SizeArgs args;

int main(int argc, char** argv) {
  if (argc != 3) {
    cerr << "Usage: <program> <step_size> <cross_tolerance>" << endl;
    return -1;
  }
  args.fixedStepSize = atof(argv[1]);
  double tol = atof(argv[2]);
  args.maxTime = 18000;
  CrossAxis<1, 4> cross;
  StandardOdeSolution sol = StandardOdeSolutionFromIvp(TwoBodies());

  auto result = ExtendPastZero(sol, RK4(), Dv(), args, cross);
  if (result != SolverResult::kOk) {
    LogResult(result);
    return -2;
  }
  cout << setprecision(10);
  auto [tperiheliom, periheliom] = CrossPoint(cross, Dv(), sol, tol, RK4());
  cout << "# periheliom\n";
  cout << "# time: " << tperiheliom << "\n";
  cout << "# point: " << periheliom[0] << " " << periheliom[1] << endl;
  result = ExtendPastZero(sol, RK4(), Dv(), args, cross);
  if (result != SolverResult::kOk) {
    LogResult(result);
    return -2;
  }
  auto [tclosing, closing] = CrossPoint(cross, Dv(), sol, tol, RK4());
  cout << "# closing point\n";
  cout << "# time: " << tclosing << "\n";
  cout << "# point: " << closing[0] << " " << closing[1] << endl;
  cout << "# initial point: " << sol.x[0][0] << " " << sol.x[0][1] << endl;
  PrintSolution(cout, sol, {0, 1}, 3000);
}

