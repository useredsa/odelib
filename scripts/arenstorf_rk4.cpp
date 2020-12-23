#include <iostream>
// Select a problem you want to solve
#include "problems/arenstorf.hpp"
// Include the method you will use in the problem
#include "methods/rk4.hpp"
#include "solvers/plain_method_solver.hpp"
// Include a container for the solution
#include "solutions/standard_ode_solution.hpp"
// Include additional tools
#include "tools/tsv_output.hpp"
using namespace std;
using namespace odelib;

SizeArgs args;

int main(int argc, char** argv) {
  if (argc != 3) {
    cerr << "Usage: <program> <max_time> <step_size>" << endl;
    return -1;
  }
  args.maxTime = atof(argv[1]);
  args.fixedStepSize = atof(argv[2]);
  StandardOdeSolution sol = StandardOdeSolutionFromIvp(Arenstorf());
  auto result = ExtendPastMaxTime(sol, RK4(), Arenstorf::Dv(), args);
  PrintSolution(cout, sol, {0, 1}, 3000);
  if (result != SolverResult::kOk) {
    return -2;
  }
}

