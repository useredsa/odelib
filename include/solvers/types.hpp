#ifndef INCLUDE_SOLVERS_HPP_
#define INCLUDE_SOLVERS_HPP_

#include <iostream>

namespace odelib {

/**
 * Uniform argument interfaces for calling the solvers
 */
struct SizeArgs {
  double tolerance;
  double fixedStepSize;
  double minStepAllowed;
  double maxStepAllowed;
  double maxTime;
};

/**
 * Enum class containing the possible outcomes
 * after calling a solver.
 */
enum class SolverResult {
  kViolatedPrecondition,
  kOk,
  kExhaustedInterval,
  kStepWentBelowMin,
  kFailedToSolveImplicitEq
};

/**
 * Logs the result of a method to the standard error output
 * Returns wheteher it was an error state.
 */
bool LogResult(SolverResult result) {
  // For when it works:
  // using enum SolverResult;
  switch (result) {
  case SolverResult::kViolatedPrecondition:
    std::cerr << "Violated precondition of method" << std::endl;
    break;
  case SolverResult::kExhaustedInterval:
    std::cerr << "Method exhausted the maximum time interval" << std::endl;
    break;
  case SolverResult::kStepWentBelowMin:
    std::cerr << "Step size went below minimum" << std::endl;
    break;
  case SolverResult::kFailedToSolveImplicitEq:
    std::cerr << "Method found it impossible to solve an implicit equation"
        << std::endl;
    break;
  case SolverResult::kOk:
    std::cerr << "Method finished succesfully" << std::endl;
    return false;
  }
  return true;
}

}  // namespace odelib

#endif  // INCLUDE_SOLVERS_HPP_

