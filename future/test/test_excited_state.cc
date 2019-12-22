#include "ExcitedState.h"
#include "TestUtilities.h"

int main(){
    
    TestUtilities test;
    ExcitedState exc;

    exc.set_excitation_energy(1e6);
    exc.set_total_width(2.);
    exc.set_ground_state_width(1.);
    exc.set_two_J(1);
    exc.set_parity(true);

    test.test_equality<double, double>(exc.get_excitation_energy(), 1e6);
    test.test_equality<double, double>(exc.get_total_width(), 2.);
    test.test_equality<double, double>(exc.get_ground_state_width(), 1.);
    test.test_equality<double, double>(exc.get_two_J(), 1);
    test.test_equality<bool, bool>(exc.get_parity(), true);
}