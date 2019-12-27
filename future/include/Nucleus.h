#pragma once

#include <vector>

#include "ExcitedState.h"

using std::vector;

/// \brief Class for a nucleus
class Nucleus{
public:
    /// \brief Default constructor.
    Nucleus() = default;
    /// \brief Default destructor.
    ~Nucleus() = default;
    /// \brief Constructor for the initialization of all members.
    /// \param exc_sta Vector of ExcitedState objects.
    /// \param m Nuclear mass.
    /// \param tJ Ground-state angular momentum times 2.
    Nucleus(vector<ExcitedState> exc_sta, const double m, const unsigned int tJ);

    /// \brief Function to add an excited state.
    /// \param exc ExcitedState object.
    void add_excited_state( ExcitedState &exc ){ excited_states.push_back(exc); };
    /// \brief Function to get the number of excited states.
    size_t get_n_excited_states() const { return excited_states.size(); };
    /// \brief Function which returns an excited state of the nucleus.
    /// \param i Number of the excited state to return.
    ExcitedState& get_excited_state( const size_t i ) { return excited_states.at(i); };

    /// \brief Function which returns the ground-state angular momentum.
    unsigned int get_two_J() const { return two_J; };
    /// \brief Function to set the ground-state angular momentum.
    /// \param tJ New value for the angular momentum.
    void set_two_J( const unsigned int tJ ){ two_J = tJ; };
    
    /// \brief Function which returns the nuclear mass.
    unsigned int get_mass() const { return mass; };
    /// \brief Function to set the nuclear mass.
    /// \param m New value for the nuclear mass.
    void set_mass( const unsigned int m ){ mass = m; };

    /// \brief Function to calculate the Doppler width of an excited state.
    /// \param i Number of the excited state.
    /// \param temperature (Effective) temperature of the material.
    ///
    /// The Doppler width is given by [Eq. (7) in \cite Metzger1959 or Eq. (4.23) in \cite Romig2015]:
    /// \f[
    ///     \Delta = \sqrt{\frac{2 k_B T}{M c^2}} E_j,
    /// \f]
    /// where \f$k_B\f$ is the Boltzmann constant, \f$T\f$ is the (effective) temperature of the material,
    /// \f$M\f$ is the nuclear mass, \f$c\f$ is the speed of light, and \f$E_j\f$ is the excitation energy of the state.
    /// In the case that the Doppler width is much larger than the width of an excited state, i.e.:
    /// \f[
    ///     \Delta \gg \Gamma,
    /// \f]
    /// the energy dependence of the absorption cross section will have the shape of a normal distribution with
    /// the standard deviation [Eq. (6) in \cite Metzger1959 or Eq. (4.26) in \cite Romig2015]:
    /// \f[
    ///     \sigma = \frac{\Delta }{\sqrt{2}}.
    /// \f]
    double doppler_width( const size_t i, const double temperature ) const;
    /// \brief Function to calculate the Doppler width for all excited states.
    /// \param temperature (Effective) temperature of the material.
    /// \return Vector of Doppler widths.
    vector<double> doppler_width( const double temperature ) const;

    /// \brief Function to calculate the energy-integrated excitation cross section for an excited state.
    /// \param i Number of the excited state.
    ///
    /// The energy-integrated cross section is given by [Eq. (12) in \cite Metzger1959 or Eq. (4.28) in \cite Romig2015]:
    /// \f[
    ///     \int_0^{\infty} \sigma_{0 \to j}(E) \mathrm{d} E = \pi^2 \left( \frac{\hbar c}{E_j} \right)^2 \frac{2J_j + 1}{2J_0 + 1} \Gamma_{0 \to j}.
    /// \f]
    /// As can be seen on the left-hand side of the equation above, it is defined as the energy integral over the
    /// cross section \f$\sigma_{0 \to j}\f$ for the absorption of a gamma ray by a nucleus in the ground state \f$0\f$ and
    /// the subsequent excitation of the nucleus to a state \f$j\f$.
    /// On the right-hand side of the equation, \f$\hbar\f$ is the reduced Planck constant, \f$c\f$ is the speed of light, \f$E_j\f$ is the excitation
    /// energy of the excited state and \f$J_j\f$ its angular momentum, \f$J_0\f$ is the ground-state angular momentum,
    /// and \f$\Gamma_{0 \to j}\f$ is the partial width for the excitation of state \f$j\f$ from the ground state.
    double energy_integrated_cs( const size_t i ) const;
    /// \brief Function to calculate the energy-integrated excitation cross section for all excited states.
    vector<double> energy_integrated_cs() const;

    /// \brief Function to create an energy grid for the evaluation of the excitation cross sections of all excited states.
    /// \param e_min Start energy of the grid.
    /// \param e_max End energy of the grid.
    /// \param n_energies_per_state Number of sampled points per excited state. In total, the grid will consists of \f$n_s n_e\f$ points, where \f$n_s\f$ is the number of excited states, and \f$n_e\f$ the number of energies per state.
    /// This function assumes that the nucleus is at rest, i.e. the cross sections will all have a Breit-Wigner form.
    vector<double> energies(const double e_min, const double e_max, const size_t n_energies_per_state) const;
    /// \brief Function to calculate the coverage of every single cross section.
    /// \param e_min Start energy of the grid.
    /// \param e_max End energy of the grid.
    vector<double> cross_section_coverage(const double e_min, const double e_max) const;
    /// \brief Function to calculate the total absorption cross section of a nucleus at the given energies.
    /// \param energies Vector of energies, at which the cross section should be evaluated.
    ///
    /// The energy-dependent total cross section for excitation of a nucleus from the ground state \f$\sigma_0(E)\f$ is given by a sum over
    /// all cross sections for the excitation from the ground state \f$0\f$ to states \f$j\f$:
    /// \f[
    ///     \sigma_0 (E) = \sum_j \sigma_{0 \to j} (E).
    /// \f]
    /// The latter have the shape of a Breit-Wigner distribution, and they are given by [Eq. (3) in \cite Metzger1959 or Eq. (4.17) in \cite Romig2015]:
    /// \f[
    ///     \sigma_{0 \to j} (E) = \frac{\pi}{2} \left( \frac{\hbar c}{E_j} \right)^2 \frac{2J_j + 1}{2J_0 + 1} \frac{\Gamma_{0 \to j} \Gamma_j}{(E - E_j)^2 + \Gamma_{j}^2/4}
    /// \f]
    /// In this equation, \f$\hbar\f$ is the reduced Planck constant, \f$c\f$ is the speed of light, \f$E_j\f$ is the excitation
    /// energy of the excited state and \f$J_j\f$ its angular momentum, \f$J_0\f$ is the ground-state angular momentum,
    /// \f$\Gamma_{0 \to j}\f$ is the partial width for the excitation of state \f$j\f$ from the ground state, and \f$\Gamma_j\f$ is the total
    /// width of the excited state.
    vector<double> cross_section(const vector<double> &energies) const;
    /// \brief Function to calculate the Doppler-broadened total absorption cross section of a nucleus at the given energies.
    /// \param energies Vector of energies at which the cross section should be evaluated.
    /// \param temperature (Effective) temperature of the material.
    ///
    /// It is assumed that, at a finite temperature \f$T\f$, the probability distribution \f$p\f$ for the velocity component of the nuclei parallel to 
    /// the incoming gamma ray \f$v_\parallel\f$, is given by a Maxwell-Boltzmann distribution [Eq. (5) in \cite Metzger1959 or Eq. (4.19) in \cite Romig2015]:
    /// \f[
    ///     p(v_\parallel) \mathrm{d} v_\parallel = \sqrt{\frac{M}{2\pi k_B T}} \exp \left( - \frac{M v_\parallel^2}{2 k_B T} \right)
    /// \f]
    /// In the equation above, \f$M\f$ denotes the nuclear mass, \f$k_B\f$ the Boltzmann constant, and \f$T\f$ the temperature of the material.
    /// The temperature corresponds to the actual thermodynamic temperature of the material only in the case of an ideal gas.
    /// In a realistic gaseous/liquid/solid material, the influence of the atomic binding can be taken into account at first order by replacing the
    /// temperature \f$T\f$ by an effective value (For a general discussion, see p. 58 in \cite Metzger1959. For a discussion of the effect in crystals, see
    /// the historical publication by Lamb \cite Lamb1939. For a general treatment of the atomic environment, see \cite Singwi1960).
    ///
    /// Consequently, the energy-dependent Doppler-broadened cross section \f$\sigma^D_{0 \to j}\f$ for the absorption of a gamma ray by a nucleus in the ground state \f$0\f$ 
    /// and the subsequent excitation to a state \f$j\f$ is given by a convolution of the Breit-Wigner shaped cross section at rest with the velocity distribution above,
    /// noting that the Doppler-shifted gamma-ray energy \f$E_\gamma (E_{\gamma, 0}, v_\parallel)\f$ is given by [Eq. (4) in \cite Metzger1959 or Eq. (4.20) in \cite Romig2015]:
    /// \f[
    ///     E_\gamma (E_{\gamma, 0}, v_\parallel) \approx \left( 1 + \frac{v_\parallel}{c}\right) E_{\gamma, 0}
    /// \f]
    /// Here, \f$E_{\gamma, 0}\f$ denotes the gamma-ray energy as it would appear to a nucleus at rest.
    /// The equation above is a good approximation if the velocity of the atoms is much smaller than the speed of light \f$c\f$.
    /// Using the expression for the Doppler-shifted resonance energy, the Doppler-broadened cross section is given by [Eq. (8) in \cite Metzger1959 or Eq. (4.22) in \cite Romig2015]:
    /// \f[
    ///     \sigma^D_{0 \to j} (E) = \int_{-\infty}^{\infty} \sigma_{0 \to j} \left[ E_\gamma(E, v_\parallel) \right] p(v_\parallel) \mathrm{d} v_\parallel
    /// \f]
    /// The integral can be transformed into an integral over the energy with the replacement:
    /// \f[
    ///     v_\parallel = \frac{c}{E_{\gamma, 0}} \left( E_\gamma - E_{\gamma, 0} \right)
    /// \f]
    /// \f[
    ///     \mathrm{d} v_\parallel = \frac{c}{E_{\gamma, 0}} \mathrm{d} E_\gamma.
    /// \f]
    /// Since \f$v_\parallel\f$ depends on the difference between the Doppler-shifted gamma-ray energy and the one at rest,
    /// this yields a convolution of a Breit-Wigner and a normal distribution:
    /// \f[
    ///     \sigma^D_{0 \to j} (E) = \int_0^\infty \sigma_{0 \to j} \left( E_\gamma \right) p \left( E_\gamma - E \right) \mathrm{d} E_\gamma.
    /// \f]
    /// The result of this convolution is a Voigt profile \cite VoigtProfile2019 \f$V(E, \sigma, \Gamma)\f$ with a standard deviation
    /// \f[
    ///     \sigma = \frac{\Delta}{\sqrt{2}},
    /// \f]
    /// and a width 
    /// \f[
    ///     \Gamma = \Gamma_j
    /// \f]
    /// In the equations above, \f$\Delta\f$ denotes the Doppler width of the excited state,
    /// and \f$\Gamma_j\f$ its total width.
    vector<double> cs_doppler_broadened(const vector<double> &energies, const double temperature) const;

protected:
    vector<ExcitedState> excited_states; ///< Vector of excited states of the nucleus.

    double mass; ///< Nuclear mass (atomic mass units times the speed of light squared, i.e. \f$uc^2\f$).
    unsigned int two_J; ///< Total angular momentum of the ground state times 2, in units of the reduced Planck constant.
};