#pragma once

#include <vector>

#include "ExcitedState.h"

using std::vector;

class Nucleus{
public:
    Nucleus();

    void add_excited_state( ExcitedState &exc ){ excited_states.push_back(exc); };
    size_t get_n_excited_states() const { return excited_states.size(); };
    ExcitedState get_excited_state( size_t i ) const { return excited_states.at(i); };

    unsigned int get_two_J() const { return two_J; };
    void set_two_J( unsigned int tJ ){ two_J = tJ; };

    vector<double> energies(const double e_min, const double e_max, const size_t n_energies) const;

private:
    vector<ExcitedState> excited_states;

    double mass; // atomic mass units u
    unsigned int two_J;
};