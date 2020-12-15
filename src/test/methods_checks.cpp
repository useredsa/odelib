#include "methods/interfaces/plain_method.hpp"
#include "methods/interfaces/adaptive_single_step_method.hpp"
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
static_assert(PlainMethod<Euler<Arenstorf::Dv>>);
static_assert(PlainMethod<ModEuler<Arenstorf::Dv>>);
static_assert(PlainMethod<RK4<Arenstorf::Dv>>);
static_assert(PlainMethod<Taylor<Taylor1::Dv, 3>>);

//TODO change name
// PlainAdaptive
static_assert(AdaptiveSingleStepMethod<
                RichardsonExtrapolation<Euler<Arenstorf::Dv>, Arenstorf::Dv>>);
static_assert(AdaptiveSingleStepMethod<Fehlberg<Arenstorf::Dv>>);

// PlainMultistep
static_assert(PlainMultistepMethod<AdamsBashforth4<Arenstorf::Dv>>);

// AdaptiveMultistep
static_assert(AdaptiveMultistepMethod<PredictorCorrector4<Arenstorf::Dv>>);

// PlainImplicit
static_assert(PlainImplicitMethod<BackwardsEuler<Arenstorf::Dv>>);
static_assert(PlainImplicitMethod<Trapezoidal<Arenstorf::Dv>>);

// Bdf
// static_assert(BackwardDifferentiationFormula<Bdf1<Rigid1::Dv>>);
static_assert(BackwardDifferentiationFormula<Bdf2<Rigid1::Dv>>);
static_assert(BackwardDifferentiationFormula<Bdf3<Rigid1::Dv>>);
static_assert(BackwardDifferentiationFormula<Bdf4<Rigid1::Dv>>);
static_assert(BackwardDifferentiationFormula<Bdf5<Rigid1::Dv>>);
static_assert(BackwardDifferentiationFormula<Bdf6<Rigid1::Dv>>);

}
