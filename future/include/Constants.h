#pragma once

/// \brief Struct to store mathematical and physical constants
///
/// Mathematical constants were evaluated using NumPy \cite Oliphant2006 \cite VanDerWalt2011.
/// Many of them could also be provided by the `<numbers>` header of C++20 \cite Numbers2019 in the future.
///
/// For each physical constant, the respective source is given.
struct Constants{
    static constexpr double pi = 3.141592653589793; ///< pi
    static constexpr double pi_squared = 9.869604401089358; ///< pi*pi
    static constexpr double inverse_pi = 0.3183098861837907; ///< 1/pi
    static constexpr double inverse_twopi = 0.15915494309189535; ///< 1/(2*pi)
    static constexpr double inverse_sqrt_two = 0.7071067811865475; ///< 1/sqrt(2)
    static constexpr double inverse_sqrt_twopi = 0.3989422804014327; ///< 1/sqrt(2*pi)

    // Physical constants
    static constexpr double kB = 8.617333262145e-5; ///< Boltzmann constant (eV K^-1) \cite BIPM2019
    static constexpr double u = 931.49410242e6; ///< Atomic mass unit (eV c^-2) \cite Dalton2019
    static constexpr double hbarc = 197.3269804e6; ///< Reduced Planck constant times the speed of light (eV fm) \cite PlanckConstant2019
    static constexpr double hbarc_squared = 3.893793719378198e+16; ///< Reduced Planck constant times the speed of light, squared (eV^2 fm^2) \cite PlanckConstant2019
};