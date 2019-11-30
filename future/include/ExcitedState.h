#pragma once

class ExcitedState{

public:
    ExcitedState(){};
    ExcitedState(const double exc_ene, const double gamma,
        const double gamma_0, const unsigned int twoJ,
        const bool p);

    double get_excitation_energy() const { return excitation_energy; };
    void set_excitation_energy(const double exc){ excitation_energy = exc; };

    double get_ground_state_width() const { return ground_state_width; };   
    void set_ground_state_width(const double g0){ ground_state_width = g0; };

    double get_total_width() const { return total_width; };
    void set_total_width(const double g){ total_width = g; };

    unsigned int get_two_J() const { return two_J; };
    void set_two_J(const unsigned int tJ){ two_J = tJ; };

    bool get_parity() const { return parity; };
    void set_parity(const bool p) { parity = p; };

private:
    double excitation_energy;
    double total_width;
    double ground_state_width;

    unsigned int two_J;
    bool parity;
};