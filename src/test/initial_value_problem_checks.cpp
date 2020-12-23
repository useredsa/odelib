#include "initial_value_problem.hpp"

#include "problems/arenstorf.hpp"
#include "problems/two_bodies.hpp"
#include "problems/rigid1.hpp"
#include "problems/taylor1.hpp"

namespace odelib {

static_assert(InitialValueProblem<Arenstorf>);
static_assert(InitialValueProblem<TwoBodies>);
static_assert(InitialValueProblem<Rigid1>);
static_assert(InitialValueProblem<Taylor1>);

static_assert(SpaceDerivableIvpDerivative<Rigid1::Dv>);
static_assert(NDerivableIvpDerivative<Taylor1::Dv>);

}  // namespace odelib
