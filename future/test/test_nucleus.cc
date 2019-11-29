#include <numeric>
#include <vector>

#include "Constants.h"
#include "ExcitedState.h"
#include "Integrator.h"
#include "Nucleus.h"
#include "TestUtilities.h"

using std::vector;

int main(){

    TestUtilities test;

    Nucleus nuc;
    nuc.set_two_J(1);
    nuc.set_mass(1.);

    double exc1_excitation_energy = 1e6;
    double exc1_ground_state_width = 1.;
    double exc1_total_width = 1.;
    unsigned int exc1_two_J = 1;

    double exc2_excitation_energy = 2e6;
    double exc2_ground_state_width = 2.;
    double exc2_total_width = 2.;
    unsigned int exc2_two_J = 2;

    ExcitedState exc1, exc2;
    exc1.set_excitation_energy(exc1_excitation_energy);
    exc1.set_ground_state_width(exc1_ground_state_width);
    exc1.set_total_width(exc1_total_width);
    exc1.set_two_J(exc1_two_J);

    exc2.set_excitation_energy(exc2_excitation_energy);
    exc2.set_ground_state_width(exc2_ground_state_width);
    exc2.set_total_width(exc2_total_width);
    exc2.set_two_J(exc2_two_J);

    nuc.add_excited_state(exc1);
    nuc.add_excited_state(exc2);

    // Test the analytical calculation of the energy-integrated cross section, which is an
    // important benchmark for numerical integrations of cross sections.

    test.is_close_relative(nuc.energy_integrated_cs(0),
        Constants::pi_squared*Constants::hbarc_squared
            /(exc1_excitation_energy*exc1_excitation_energy), test.num_tol_rel);

    // Test calculation of Doppler width
    // For a simple result, enter a temperature for which the square-root term equals 1
    vector<double> doppler_widths = nuc.doppler_width( 0.5*nuc.get_mass()*Constants::u / Constants::kB );
    test.is_close_relative(doppler_widths[0], nuc.get_excited_state(0).get_excitation_energy(), 
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
    nuc.set_two_J(exc1.get_two_J());
    nuc.add_excited_state(exc1);
    double e_min = exc1_excitation_energy - 5*exc1_total_width;
    double e_max = exc1_excitation_energy + 5*exc1_total_width;
    vector<double> energies = nuc.energies(e_min, e_max, 2000); // Empirically determined this number of 
        //points to achieve the required numerical precision
    vector<double> cs_cov = nuc.cross_section_coverage(e_min, e_max);
    vector<double> cs = nuc.cross_section(energies);
    vector<double> cs_ene_int = nuc.energy_integrated_cs();
    test.is_close_relative(inte.trapezoidal_rule(energies, cs)/cs_cov[0],
        std::accumulate(cs_ene_int.begin(), cs_ene_int.end(), 0.), test.num_tol_rel);
}