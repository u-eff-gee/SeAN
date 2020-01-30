#include "ExcitedState.h"

ExcitedState::ExcitedState(const double exc_ene, const double gamma,
    const double gamma_0, const unsigned int tJ, const bool p):
    excitation_energy(exc_ene),
    total_width(gamma),
    ground_state_width(gamma_0),
    parity(p), two_J(tJ)
    {};