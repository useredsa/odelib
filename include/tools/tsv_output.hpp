#ifndef INCLUDE_TOOLS_TSV_OUTPUT_HPP_
#define INCLUDE_TOOLS_TSV_OUTPUT_HPP_

#include <cassert>
#include <iostream>
#include <vector>
#include "ode_solution.hpp"

namespace odelib {

/**
 * Prints the desired indices of a solution with tsv format.
 * It can also select a maximum number of points,
 * choosing only a handful of equispaced points.
 */
template <OdeSolution Sol>
void PrintSolution(std::ostream& out, const Sol& sol,
    const std::vector<int>& indices, int maxNum = 1e4) {
  assert(!indices.empty());
  if (sol.empty()) return;
  for (int id : indices) {
    assert(0 <= id && id < (int) sol.x[0].size());
  }

  long long p1 = 0, p2 = 0;
  for (size_t i = 0; i < sol.size(); ++i) {
    if (p1 >= p2) {
      out << sol.x[i][indices[0]];
      for (size_t j = 1; j < indices.size(); ++j) {
        out << '\t' << sol.x[i][indices[j]];
      }
      out << '\n';
      p2 += sol.size();
    }
    p1 += maxNum;
  }
}

/**
 * Prints the time and the desired indices of a solution with tsv format.
 * It can also select a maximum number of points,
 * choosing only a handful of equispaced points.
 */
template <OdeSolution Sol>
void PrintSolutionWithTime(std::ostream& out, const Sol& sol,
    const std::vector<int>& indices, int maxNum = 1e4) {
  if (sol.empty()) return;
  for (int id : indices) {
    assert(0 <= id && id < (int) sol.x[0].size());
  }

  long long p1 = 0, p2 = 0;
  for (size_t i = 0; i < sol.size(); ++i) {
    if (p1 >= p2) {
      out << sol.t[i];
      for (size_t j = 0; j < indices.size(); ++j) {
        out << '\t' << sol.x[i][indices[j]];
      }
      out << '\n';
      p2 += sol.size();
    }
    p1 += maxNum;
  }
}

template <OdeSolution Sol>
void PrintSolution(std::ostream& out, const Sol& sol, int maxNum = 1e4) {
  if (sol.empty()) return;
  std::vector<int> ind(sol.x[0].size());
  for (int i = 0; i < sol.x[0].size(); ++i) {
    ind[i] = i;
  }
  PrintSolution(out, sol, ind, maxNum);
}

template <OdeSolution Sol>
void PrintSolutionWithTime(std::ostream& out, const Sol& sol, int maxNum = 1e4) {
  if (sol.empty()) return;
  std::vector<int> ind(sol.x[0].size());
  for (int i = 0; i < sol.x[0].size(); ++i) {
    ind[i] = i;
  }
  PrintSolutionWithTime(out, sol, ind, maxNum);
}

}  // namespace odelib

#endif  // INCLUDE_TOOLS_TSV_OUTPUT_HPP_

