#include "Constants.h"
#include "Target.h"
#include "TestUtilities.h"

int main(){
    TestUtilities test;

    ExcitedState exc1;
    exc1.set_excitation_energy(1e6);
    exc1.set_ground_state_width(1.);
    exc1.set_total_width(1.);
    exc1.set_two_J(1);
    ExcitedState exc2;
    exc2.set_excitation_energy(1e6+10.);
    exc2.set_ground_state_width(1.);
    exc2.set_total_width(1.);
    exc2.set_two_J(1);

    Nucleus nuc;
    nuc.set_mass(1.);
    nuc.set_two_J(1);
    nuc.add_excited_state(exc1);
    nuc.add_excited_state(exc2);

    Target target;
    target.set_temperature(100.);
    
}