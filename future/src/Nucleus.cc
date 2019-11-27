#include <algorithm>

#include "Math/QuantFuncMathCore.h"
#include "Math/ProbFuncMathCore.h"

#include "Nucleus.h"

Nucleus::Nucleus():excited_states(0){}

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