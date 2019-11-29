#pragma once

#include "Nucleus.h"

class Target{

public:
    Target() = default;

    void add_nucleus( Nucleus &exc ){ nuclei.push_back(exc); };
    size_t get_n_nuclei() const { return nuclei.size(); };
    Nucleus get_nucleus( size_t i ) const { return nuclei.at(i); };

    unsigned int get_temperature() const { return temperature; };
    void set_temperature( unsigned int t ){ temperature = t; };

    vector<double> energies(const double e_min, const double e_max, const size_t n_energies) const;

private:
    vector<Nucleus> nuclei;

    double temperature; // K
};