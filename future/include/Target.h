#pragma once

#include "Nucleus.h"

/// \brief Class for a target
///
/// A 'target' means a piece of material which is hit by a gamma-ray beam.
/// It may consist of several different species of nuclei.
class Target{

public:
    /// \brief Default constructor.
    Target() = default;
    /// \brief Default destructor.
    ~Target() = default;
    /// \brief Constructor for the initialization of all members.
    /// \param nuc Vector of Nucleus objects.
    /// \param tem Temperature of the target.
    Target(vector<Nucleus> nuc, const double tem);

    /// \brief Function to add a nucleus.
    /// \param nuc Nucleus object.
    void add_nucleus( Nucleus &nuc ){ nuclei.push_back(nuc); };
    /// \brief Function to get the number of nuclei.
    size_t get_n_nuclei() const { return nuclei.size(); };
    /// \brief Function which returns a nucleus of the target.
    /// \param i Number of the nucleus to return.
    Nucleus& get_nucleus( size_t i ) { return nuclei.at(i); };

    /// \brief Function to get the temperature of the target.
    double get_temperature() const { return temperature; };
    /// \brief Function to set the temperature of the target.
    /// \param t New value of the temperature.
    void set_temperature( const double t ){ temperature = t; };

    /// \brief Function to create an energy grid for the evaluation of the excitation cross sections of all excited states of all nuclei.
    /// \param e_min Start energy of the grid.
    /// \param e_max End energy of the grid.
    /// \param n_energies_per_state Number of sampled points per excited state.
    /// In total, the grid will consist of 
    /// \f[
    ///     \sum_n n_{n,s} n_e
    /// \f] 
    /// points, where the sum is assumed to run over all nuclei \f$n\f$. 
    /// The symbol \f$n_{n,s}\f$ denotes the number of excited states of the nucleus \f$n\f$, 
    /// and \f$n_e\f$ is the number of energies per state.
    vector<double> energies(const double e_min, const double e_max, const size_t n_energies_per_state);

    /// \brief Function to call the total Doppler-broadened absorption cross section of the target.
    /// \param energies Vector of energies at which the cross section should be evaluated.
    vector<double> cross_section(const vector<double> &energies) const;

protected:
    vector<Nucleus> nuclei; ///< Vector of nuclei of the target.

    double temperature; ///< Temperature of the target (K).
};