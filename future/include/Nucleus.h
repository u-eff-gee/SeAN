#pragma once

#include <vector>

#include "ExcitedState.h"

using std::vector;

class Nucleus{
public:
    Nucleus():excited_states(0){};

    void add_excited_state( ExcitedState &exc ){ excited_states.push_back(exc); };
    size_t get_n_excited_states() const { return excited_states.size(); };
    ExcitedState get_excited_state(size_t i) const { return excited_states.at(i); };

private:
    vector<ExcitedState> excited_states;
    double mass;
    unsigned int two_J;
};