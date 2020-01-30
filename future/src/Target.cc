#include "Constants.h"
#include "Grid.h"
#include "Target.h"

Target::Target(vector<Nucleus> nuc, const double tem):
    nuclei(nuc), temperature(tem)
    {}

vector<double> Target::energies(const double e_min, const double e_max,
    const size_t n_energies_per_state) {
    size_t n_energies = 0;
    for(size_t i = 0; i < nuclei.size(); ++i){
        n_energies += n_energies_per_state*nuclei[i].get_n_excited_states();
    }

    vector<double> ene(n_energies, 0.);

    vector<double> ene_state;
    Grid grid;
    unsigned int index = 0;

    for(size_t i = 0; i < nuclei.size(); ++i){
        for(size_t j = 0; j < nuclei[i].get_n_excited_states(); ++j){
            if(nuclei[i].doppler_width(i, temperature) > nuclei[i].get_excited_state(j).get_total_width()){
                ene_state = grid.normal(nuclei[i].get_excited_state(j).get_excitation_energy(),
                    nuclei[i].doppler_width(j, temperature)*Constants::inverse_sqrt_two,
                    e_min, e_max, n_energies_per_state);
            } else{
                ene_state = grid.breit_wigner(nuclei[i].get_excited_state(j).get_excitation_energy(),
                    nuclei[i].get_excited_state(j).get_total_width(),
                    e_min, e_max, n_energies_per_state);                
            }
            for(size_t k = 0; k < n_energies_per_state; ++k){
                ene[k+index] = ene_state[k];
            }
            index += n_energies_per_state;
        }
    }
    grid.strictly_increasing(ene);
    return ene;
}

vector<double> Target::cross_section(const vector<double> &energies) const {
    vector<double> cs(energies.size(), 0.);
    vector<double> cs_nuc(energies.size());

    for(size_t i = 0; i < nuclei.size(); ++i){
        cs_nuc = nuclei[i].cs_doppler_broadened(energies, temperature);
        for(size_t j = 0; j < energies.size(); ++j){
            cs[j] += cs_nuc[j];
        }
    }

    return cs;
}