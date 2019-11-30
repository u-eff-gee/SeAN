#include "Constants.h"
#include "ExcitedState.h"
#include "Integrator.h"
#include "Target.h"
#include "TestUtilities.h"

int main(){
    TestUtilities test;

    vector<ExcitedState> exc{ExcitedState(1e6-5., 1., 1., 1, 1), 
        ExcitedState(1e6+5., 1., 1., 1, 1)};
    
    // With the given choices of the nuclear mass and the temperature,
    // the Doppler width will be equal to 1e-6 times the excitation energy
    vector<Nucleus> nuc{Nucleus(exc, 1./Constants::u, 1)};
    double temperature = 1e-12/(2.*Constants::kB);
    test.is_equal<double, double>(nuc[0].doppler_width(0, temperature), (1e6-5.)*1e-6);
    test.is_equal<double, double>(nuc[0].doppler_width(1, temperature), (1e6+5.)*1e-6);

    // Set the total widths of the first (second) excited state smaller (larger) than
    // the Doppler width to obtain a more Breit-Wigner (normal-) like distribution.
    nuc[0].get_excited_state(0).set_total_width(2.0*nuc[0].doppler_width(0, temperature));
    nuc[0].get_excited_state(0).set_ground_state_width(nuc[0].get_excited_state(0).get_total_width());
    nuc[0].get_excited_state(1).set_total_width(0.5*nuc[0].doppler_width(1, temperature));
    nuc[0].get_excited_state(1).set_ground_state_width(nuc[0].get_excited_state(1).get_total_width());

    Target target(nuc, temperature);
    test.is_equal<double, double>(temperature, target.get_temperature());

    test.larger_than(target.get_nucleus(0).get_excited_state(0).get_total_width(),
        target.get_nucleus(0).doppler_width(0, target.get_temperature()));
    test.larger_than(target.get_nucleus(0).doppler_width(1, target.get_temperature()),
        target.get_nucleus(0).get_excited_state(1).get_total_width());

    double energy_range = 500.;
    vector<double> energies = target.energies(1e6-energy_range, 1e6+energy_range, 5);
    test.is_equal<size_t, size_t>(energies.size(), 10);

    // Use a more strict numerical tolerance here, because differences on the order of
    // 1 eV are compared to energies of 1 MeV
    test.is_equal<double, double>(energies[0], 1e6-energy_range);

    test.is_close_relative(energies[1],
        target.get_nucleus(0).get_excited_state(0).get_excitation_energy()-
        0.5*target.get_nucleus(0).get_excited_state(0).get_total_width()-
        0.5*(target.get_nucleus(0).get_excited_state(0).get_excitation_energy()-
        0.5*target.get_nucleus(0).get_excited_state(0).get_total_width()-1e6+energy_range), 1e-2*test.num_tol_rel);
    test.is_close_relative(energies[2],
        target.get_nucleus(0).get_excited_state(0).get_excitation_energy()
        -0.5*target.get_nucleus(0).get_excited_state(0).get_total_width(),
        1e-2*test.num_tol_rel);
    test.is_close_relative(energies[3],
        target.get_nucleus(0).get_excited_state(0).get_excitation_energy(), 1e-2*test.num_tol_rel);
    test.is_close_relative(energies[4],
        target.get_nucleus(0).get_excited_state(0).get_excitation_energy()
        +0.5*target.get_nucleus(0).get_excited_state(0).get_total_width(),
        1e-2*test.num_tol_rel);

    test.is_close_relative(energies[5],
        target.get_nucleus(0).get_excited_state(1).get_excitation_energy()-
        0.6744897501960817*Constants::inverse_sqrt_two
        *target.get_nucleus(0).doppler_width(1, temperature), 1e-2*test.num_tol_rel);
    test.is_close_relative(energies[6],
        target.get_nucleus(0).get_excited_state(1).get_excitation_energy(), 1e-2*test.num_tol_rel);
    test.is_close_relative(energies[7],
        target.get_nucleus(0).get_excited_state(1).get_excitation_energy()+
        0.6744897501960817*Constants::inverse_sqrt_two
        *target.get_nucleus(0).doppler_width(1, temperature), 1e-2*test.num_tol_rel);
    test.is_close_relative(energies[8],
        1e6+energy_range - 0.5*(1e6+energy_range-
        (target.get_nucleus(0).get_excited_state(1).get_excitation_energy()+
        0.6744897501960817*Constants::inverse_sqrt_two
        *target.get_nucleus(0).doppler_width(1, temperature))), 1e-2*test.num_tol_rel);

    test.is_equal<double, double>(energies[9], 1e6+energy_range);

    // Test the calculation of the Doppler-broadened resonance by comparing the
    // numerical integral to the analytical solution.
    // This is done for several values of the temperature to study the transition
    // from a Breit-Wigner to a normal distribution.
    Integrator inte;
    target = Target(
        vector<Nucleus>{
            Nucleus(vector<ExcitedState>{
                ExcitedState(1e6, 1., 1., 1, true)
            }, 1./Constants::u, 1)
        }, 1.
    );

    vector<double> temperatures{1e-4, 1e-2, 1., 1e2};
    for(size_t i = 0; i < temperatures.size(); ++i){
        temperatures[i] = 1e-12/(2.*Constants::kB)*temperatures[i];
    } // Corresponds to Doppler widths 0.01, 0.1, 1., and 10.

    double coverage = 1.;
    //for(auto t: temperatures){
        double t = temperatures[3];
        coverage = 1.;

        target.set_temperature(t);
        double range = 20.;
        energies = target.energies(1e6-range, 1e6+range, 10000);
        vector<double> cs = target.cross_section(energies);
        if(target.get_nucleus(0).doppler_width(0, t) <= target.get_nucleus(0).get_excited_state(0).get_total_width())
            coverage = target.get_nucleus(0).cross_section_coverage(energies[0], energies[energies.size()-1])[0];
        test.is_close_relative(inte.trapezoidal_rule(energies, cs),
            target.get_nucleus(0).energy_integrated_cs(0)*coverage,
            test.num_tol_rel);
    //}
 }