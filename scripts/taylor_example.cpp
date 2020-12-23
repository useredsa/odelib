#include <iostream>
// Select a problem you want to solve
#include "problems/taylor1.hpp"
// Include the method you will use in the problem
#include "methods/taylor.hpp"
#include "solvers/plain_method_solver.hpp"
// Include a container for the solution
#include "solutions/standard_ode_solution.hpp"
// Include additional tools
#include "tools/tsv_output.hpp"
using namespace std;
using namespace odelib;

SizeArgs args;

template <int Order>
SolverResult solve(StandardOdeSolution<1>& sol) {
  return ExtendPastMaxTime(sol, Taylor<Order>(), Taylor1::Dv(), args);
}

int main(int argc, char** argv) {
  if (argc != 4) {
    cerr << "Usage: <program> <max_time> <step_size> <1 <= order <= 8>" << endl;
    return -1;
  }
  args.maxTime = atof(argv[1]);
  args.fixedStepSize = atof(argv[2]);
  int order = atoi(argv[3]);
  StandardOdeSolution sol = StandardOdeSolutionFromIvp(Taylor1());

  SolverResult result = SolverResult::kOk;
  switch (order) {
    case 1:
      result = solve<1>(sol);
      break;
    case 2:
      result = solve<2>(sol);
      break;
    case 3:
      result = solve<3>(sol);
      break;
    case 4:
      result = solve<4>(sol);
      break;
    case 5:
      result = solve<5>(sol);
      break;
    case 6:
      result = solve<6>(sol);
      break;
    case 7:
      result = solve<7>(sol);
      break;
    case 8:
      result = solve<8>(sol);
      break;
  }

  double absErr = AbsDiff(sol, Taylor1::AnalyticalSolution());
  double meanErr = MeanDiff(sol, Taylor1::AnalyticalSolution());
  cout << "# Absolute Error: " << absErr << "\n";
  cout << "# Mean Error: " << meanErr << "\n";
  PrintSolution(cout, sol, 3000);
  if (result != SolverResult::kOk) {
    return -2;
  }
}

