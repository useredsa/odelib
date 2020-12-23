#include <iostream>
// Select a problem you want to solve
#include "problems/arenstorf.hpp"
// Include the method you will use in the problem
#include "methods/predictor_corrector_4.hpp"
#include "solvers/solve.hpp"
// Include a container for the solution
#include "solutions/standard_ode_solution.hpp"
// Include additional tools
#include "tools/tsv_output.hpp"
using namespace std;
using namespace odelib;

SizeArgs args;

int main(int argc, char** argv) {
  if (argc != 5) {
    cerr << "Usage: <program> <max_time> <min_step_allowed> <max_step_allowed> "
        << "<tolerance>" << endl;
    return -1;
  }
  args.maxTime = atof(argv[1]);
  args.minStepAllowed = atof(argv[2]);
  args.maxStepAllowed = atof(argv[3]);
  args.tolerance = atof(argv[4]);
  auto sol = SolvePastMaxTime(Arenstorf(), PredictorCorrector4(), args);
  PrintSolution(cout, sol, {0, 1}, 3000);
}

