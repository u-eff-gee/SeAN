#pragma once

/// \brief Class for an excited state of a nucleus
class ExcitedState{

public:
    /// \brief Default constructor.
    ExcitedState() = default;
    /// \brief Default destructor.
    ~ExcitedState() = default;
    /// \brief Constructor for the initialization of all members.
    /// \param exc_ene Excitation energy.
    /// \param gamma Total width.
    /// \param gamma_0 Partial decay width to the ground state.
    /// \param twoJ Total angular momentum times 2.
    /// \param p Parity.
    ExcitedState(const double exc_ene, const double gamma,
        const double gamma_0, const unsigned int twoJ,
        const bool p);

    /// \brief Function to get the excitation energy.
    double get_excitation_energy() const { return excitation_energy; };
    /// \brief Function to set the excitation energy.
    /// \param exc New value for the excitation energy.
    void set_excitation_energy(const double exc){ excitation_energy = exc; };

    /// \brief Function to get the ground-state decay width.
    double get_ground_state_width() const { return ground_state_width; };   
    /// \brief Function to set the ground-state decay width.
    /// \param g0 New value for the ground-state decay width.
    void set_ground_state_width(const double g0){ ground_state_width = g0; };

    /// \brief Function to get the total width.
    double get_total_width() const { return total_width; };
    /// \brief Function to set the total width.
    /// \param g New value for the total width.
    void set_total_width(const double g){ total_width = g; };

    /// \brief Function to get the total angular momentum.
    unsigned int get_two_J() const { return two_J; };
    /// \brief Function to set the total angular momentum.
    /// \param tJ New value for the angular momentum.
    void set_two_J(const unsigned int tJ){ two_J = tJ; };

    /// \brief Function to get the parity.
    bool get_parity() const { return parity; };
    /// \brief Function to set the parity.
    /// \param p New value for the parity.
    void set_parity(const bool p) { parity = p; };

protected:
    double excitation_energy; ///< Excitation energy (eV).
    double total_width; ///< Total width (eV).
    double ground_state_width; ///< Partial width for the decay to the ground state (eV).

    /// \brief Parity
    /// * true, 1 = positive parity
    /// * false, 0 = negative parity
    bool parity;
    unsigned int two_J; ///< Total angular momentum times 2, in units of the reduced Planck constant.
};