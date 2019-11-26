#pragma once

class ExcitedState{

public:
    ExcitedState(){};

    void set_excitation_energy(double exc){ excitation_energy = exc; };
    void set_ground_state_width(double g0){ ground_state_width = g0; };
    void set_total_width(double g){ total_width = g; };

    void set_two_J(unsigned int tJ){ two_J = tJ; };

    double get_excitation_energy() const { return excitation_energy; };
    double get_ground_state_width() const { return ground_state_width; };
    double get_total_width() const { return total_width; };

    unsigned int get_two_J() const { return two_J; };

private:
    double excitation_energy;
    double ground_state_width;
    double total_width;

    unsigned int two_J;
};