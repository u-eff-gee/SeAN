#include <numeric>
#include <vector>

#include "Constants.h"
#include "ExcitedState.h"
#include "Grid.h"
#include "Integrator.h"
#include "Nucleus.h"
#include "TestUtilities.h"

using std::vector;

int main(){

    TestUtilities test;

    vector<ExcitedState> exc{ExcitedState(1e6, 2., 1., 1, true), ExcitedState(2e6, 2., 2., 2, true)};

    Nucleus nuc{exc, 2., 2};
    nuc.set_two_J(1);
    nuc.set_mass(1.);
    nuc.get_excited_state(0).set_total_width(1.);
    test.test_equality<double, double> (nuc.get_excited_state(0).get_total_width(), 1.);

    test.test_equality<double, double>(nuc.get_mass(), 1.);
    test.test_equality<unsigned int, unsigned int>(nuc.get_two_J(), 1);

    // Test the analytical calculation of the energy-integrated cross section, which is an
    // important benchmark for numerical integrations of cross sections.

    test.test_closeness(nuc.energy_integrated_cs(0),
        Constants::pi_squared*Constants::hbarc_squared
            /(exc[0].get_excitation_energy()*exc[0].get_excitation_energy()), test.num_tol_rel);

    // Test calculation of Doppler width
    // For a simple result, enter a temperature for which the square-root term equals 1
    vector<double> doppler_widths = nuc.doppler_width( 0.5*nuc.get_mass()*Constants::u / Constants::kB );
    test.test_closeness(doppler_widths[0], nuc.get_excited_state(0).get_excitation_energy(), 
        test.num_tol_rel);

    // Calculate the cross section and test whether it is correct by
    // comparing the numerical integral to the analytical value.
    //
    // To compensate for a too-small integration range of the very slowly
    // decaying BW distribution, scale the numerical result by the inverse
    // of the expected quantile.
    // The remaining deviation will then be purely numerical.
    Integrator inte;
    nuc = Nucleus();
    nuc.set_two_J(exc[0].get_two_J());
    nuc.add_excited_state(exc[0]);
    double e_min = exc[0].get_excitation_energy() - 5*exc[0].get_total_width();
    double e_max = exc[0].get_excitation_energy() + 5*exc[0].get_total_width();
    vector<double> energies = nuc.energies(e_min, e_max, 2000); // Empirically determined this number of 
        //points to achieve the required numerical precision
    vector<double> cs_cov = nuc.cross_section_coverage(e_min, e_max);
    vector<double> cs = nuc.cross_section(energies);
    vector<double> cs_ene_int = nuc.energy_integrated_cs();
    test.test_closeness(inte.trapezoidal_rule(energies, cs)/cs_cov[0],
        std::accumulate(cs_ene_int.begin(), cs_ene_int.end(), 0.), test.num_tol_rel);
}