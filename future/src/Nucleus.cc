#include <algorithm>

#include "Math/PdfFuncMathCore.h"
#include "Math/ProbFuncMathCore.h"
#include "Math/QuantFuncMathCore.h"

#include "Constants.h"
#include "Nucleus.h"

Nucleus::Nucleus():excited_states(0){}

double Nucleus::energy_integrated_cs( const size_t i ) const {
    return Constants::pi_squared*Constants::hbarc_squared
        /(excited_states[i].get_excitation_energy()*excited_states[i].get_excitation_energy())
        *(excited_states[i].get_two_J()+1)/(two_J+1)*excited_states[i].get_ground_state_width();
}

vector<double> Nucleus::energy_integrated_cs() const {
    vector<double> ene_int_cs(excited_states.size());

    for(size_t i = 0; i < excited_states.size(); ++i){
        ene_int_cs[i] = energy_integrated_cs(i);
    }

    return ene_int_cs;
}

vector<double> Nucleus::energies(const double e_min, const double e_max,
const size_t n_energies_per_state) const 
{
    size_t n_energies = n_energies_per_state*excited_states.size()+2;
    vector<double> ene(n_energies, 0.);
    ene[0] = e_min;
    ene[n_energies-1] = e_max;

    double quantile_increment{0.}, quantile_min{0.}, quantile_max{0.};

    for(size_t i = 0; i < excited_states.size(); ++i){

        quantile_min = ROOT::Math::breitwigner_cdf(
            e_min, excited_states[i].get_total_width(), excited_states[i].get_excitation_energy());
        quantile_max = ROOT::Math::breitwigner_cdf(
            e_max, excited_states[i].get_total_width(), excited_states[i].get_excitation_energy());
        quantile_increment = (quantile_max-quantile_min)/(n_energies_per_state+1);

        for(size_t j = 1; j < n_energies_per_state + 1; ++j){
            ene[i*n_energies_per_state+j] = ROOT::Math::breitwigner_quantile(
                quantile_min+j*quantile_increment, excited_states[i].get_total_width())
                +excited_states[i].get_excitation_energy();
        }
    }

    std::sort(ene.begin(), ene.end());

    return ene;
}

vector<double> Nucleus::cross_section_coverage(const double e_min, const double e_max) const {
    vector<double> cs_cov(excited_states.size());

    for(size_t i = 0; i < excited_states.size(); ++i){
        cs_cov[i] = ROOT::Math::breitwigner_cdf(e_max,
            excited_states[i].get_total_width(), excited_states[i].get_excitation_energy())
            -ROOT::Math::breitwigner_cdf(e_min,
            excited_states[i].get_total_width(), excited_states[i].get_excitation_energy());
    }

    return cs_cov;
}

vector<double> Nucleus::cross_section(const vector<double> &energies) const {
    vector<double> cs(energies.size(), 0.);
    vector<double> ene_int_cs = energy_integrated_cs();

    for(size_t i = 0; i < excited_states.size(); ++i){
        for(size_t j = 0; j < energies.size(); ++j){
            cs[j] += ene_int_cs[i]*ROOT::Math::breitwigner_pdf(energies[j],
                excited_states[i].get_total_width(), excited_states[i].get_excitation_energy());
        }
    }

    return cs;
}