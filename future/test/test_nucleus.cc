#include <iostream>

#include "ExcitedState.h"
#include "Nucleus.h"
#include "TestUtilities.h"

int main(){

    TestUtilities test;

    Nucleus nuc;
    ExcitedState exc1, exc2;
    exc1.set_excitation_energy(1.);
    exc1.set_ground_state_width(1.);
    exc1.set_total_width(1.);
    exc1.set_two_J(1);

    exc2.set_excitation_energy(2.);
    exc2.set_ground_state_width(2.);
    exc2.set_total_width(2.);
    exc2.set_two_J(2);

    nuc.add_excited_state(exc1);
    nuc.add_excited_state(exc2);

    test.is_equal<size_t, size_t>(nuc.get_n_excited_states(), 2);

    test.is_equal<double, double>(nuc.get_excited_state(0).get_excitation_energy(), 1.);
    test.is_equal<double, double>(nuc.get_excited_state(0).get_ground_state_width(), 1.);
    test.is_equal<double, double>(nuc.get_excited_state(0).get_total_width(), 1.);
    test.is_equal<unsigned int, unsigned int>(nuc.get_excited_state(0).get_two_J(), 1);

    test.is_equal<double, double>(nuc.get_excited_state(1).get_excitation_energy(), 2.);
    test.is_equal<double, double>(nuc.get_excited_state(1).get_ground_state_width(), 2.);
    test.is_equal<double, double>(nuc.get_excited_state(1).get_total_width(), 2.);
    test.is_equal<unsigned int, unsigned int>(nuc.get_excited_state(1).get_two_J(), 2);
    
}