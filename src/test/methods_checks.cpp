#include "methods/interfaces/plain_method.hpp"
#include "methods/interfaces/plain_adaptive_method.hpp"
#include "methods/interfaces/plain_multistep_method.hpp"
#include "methods/interfaces/adaptive_multistep_method.hpp"
#include "methods/interfaces/plain_implicit_method.hpp"
#include "methods/interfaces/backward_differentiation_formula.hpp"

#include "methods/euler.hpp"
#include "methods/mod_euler.hpp"
#include "methods/rk4.hpp"
#include "methods/taylor.hpp"
#include "methods/richardson_extrapolation.hpp"
#include "methods/fehlberg.hpp"
#include "methods/adams_bashforth_4.hpp"
#include "methods/predictor_corrector_4.hpp"
#include "methods/backwards_euler.hpp"
#include "methods/trapezoidal.hpp"
#include "methods/backward_differentiation_formulas.hpp"

#include "problems/arenstorf.hpp"
#include "problems/taylor1.hpp"
#include "problems/rigid1.hpp"

namespace odelib {

// Plain
static_assert(PlainMethod<Euler>);
static_assert(PlainMethod<ModEuler>);
static_assert(PlainMethod<RK4>);
static_assert(PlainMethod<Taylor<3>>);

// PlainAdaptive
static_assert(PlainAdaptiveMethod<RichardsonExtrapolation<Euler>>);
static_assert(PlainAdaptiveMethod<Fehlberg>);

// PlainMultistep
static_assert(PlainMultistepMethod<AdamsBashforth4>);

// AdaptiveMultistep
//TODO I think we should be able to write it without `AdamsBashforth4`
// because it's the default template
static_assert(AdaptiveMultistepMethod<PredictorCorrector4<AdamsBashforth4>>);

// PlainImplicit
static_assert(PlainImplicitMethod<BackwardsEuler>);
static_assert(PlainImplicitMethod<Trapezoidal>);

// Bdf
// static_assert(BackwardDifferentiationFormula<Bdf1<Rigid1::Dv>>);
static_assert(BackwardDifferentiationFormula<Bdf2>);
static_assert(BackwardDifferentiationFormula<Bdf3>);
static_assert(BackwardDifferentiationFormula<Bdf4>);
static_assert(BackwardDifferentiationFormula<Bdf5>);
static_assert(BackwardDifferentiationFormula<Bdf6>);

}
