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

    double exc1_excitation_energy = 1.;
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
        Constants::pi_squared*Constants::hbarc_squared, test.num_tol_rel);

    // Test the method of Nucleus which creates an energy grid for the evaluation of multiple
    // Breit-Wigner (BW) like resonances.
    // In order to have a grid which is more dense close to the resonance and less dense
    // far away from it, the energy values are determined from the inverse CDF of the
    // BW distribution.
    //
    // Here, the limits of the energy interval are set far away from the actual resonances,
    // so that the CDF of both BW distributions will be 0 and 1 at both limits.
    // Three energy values within this interval are requested. They are the
    // 25%, 50%, and 75% quantiles of the BW distribution, which are at 
    //  excitation_energy - 0.5*total_width
    //  excitation_energy
    //  excitation_energy + 0.5*total_width
   
    vector<double> energies = nuc.energies(
        exc1.get_excitation_energy()
        -0.5*(exc2.get_excitation_energy()-exc1.get_excitation_energy()),
        exc2.get_excitation_energy()
        +0.5*(exc2.get_excitation_energy()-exc1.get_excitation_energy()), 3);

    test.is_equal<size_t, size_t>(energies.size(), 8);
    test.is_equal<double, double>(energies[0], exc1.get_excitation_energy()
            -0.5*(exc2.get_excitation_energy()-exc1.get_excitation_energy()));

    test.is_close_absolute(energies[1], exc1.get_excitation_energy()-0.5*exc1.get_total_width(), 
        test.num_tol_rel*exc1.get_total_width());
    test.is_close_absolute(energies[2], exc1.get_excitation_energy(), test.num_tol_rel*exc1.get_total_width());
    test.is_close_absolute(energies[3], exc1.get_excitation_energy()+0.5*exc1.get_total_width(), 
        test.num_tol_rel*exc1.get_total_width());

    test.is_close_absolute(energies[4], exc2.get_excitation_energy()-0.5*exc2.get_total_width(), 
        test.num_tol_rel*exc2.get_total_width());
    test.is_close_absolute(energies[5], exc2.get_excitation_energy(), test.num_tol_rel*exc2.get_total_width());
    test.is_close_absolute(energies[6], exc2.get_excitation_energy()+0.5*exc2.get_total_width(), 
        test.num_tol_rel*exc2.get_total_width());

    test.is_equal<double, double>(energies[7], exc2.get_excitation_energy()
        +0.5*(exc2.get_excitation_energy()-exc1.get_excitation_energy()));

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
    energies = nuc.energies(e_min, e_max, 2000); // Empirically determined this number of 
        //points to achieve the required numerical precision
    vector<double> cs_cov = nuc.cross_section_coverage(e_min, e_max);
    vector<double> cs = nuc.cross_section(energies);
    vector<double> cs_ene_int = nuc.energy_integrated_cs();
    test.is_close_relative(inte.trapezoidal_rule(energies, cs)/cs_cov[0],
        std::accumulate(cs_ene_int.begin(), cs_ene_int.end(), 0.), test.num_tol_rel);
}