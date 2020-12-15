#include <fstream>
#include "methods/fixed_step.hpp"
#include "methods/adaptive.hpp"
#include "methods/multistep.hpp"
#include "methods/adaptive_multistep.hpp"
#include "solvers.hpp"
#include "problems/arenstorf.hpp"
using namespace odelib;
using namespace std;

using Dv = Arenstorf::Dv;

template<AbstractNumericalSolution Ns>
void print_solution(ostream& out, const Ns& ns, int maxNum = 1e4) {
  long long p1 = 0, p2 = 0;
  for (size_t i = 0; i < ns.size(); ++i) {
    if (p1 >= p2) {
      out << ns.x()[i][0] << '\t' << ns.x()[i][1] << '\n';
      p2 += ns.size();
    }
    p1 += maxNum;
  }
}

SolverKwArgs kwargs;
double max_t;

int main(int argc, char** argv) {
  if (argc != 3) {
    cerr << "Usage: <program> <max_time> <step_size>" << endl;
  }
  max_t = atof(argv[1]);
  kwargs.fixedStepSize = atof(argv[2]);
  auto sol = solve<RK4<Dv>>(Arenstorf(), max_t, kwargs);
  print_solution(cout, sol, 3000);
}

