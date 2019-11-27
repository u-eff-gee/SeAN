#pragma once

struct Constants{
    
    // Mathematical constants (evaluated using numpy)
    // Are also provided by the <numbers> header of c++20
    static constexpr double pi = 3.141592653589793;
    static constexpr double inverse_pi = 0.3183098861837907;
    static constexpr double inverse_sqrt_twopi = 0.3989422804014327;

    // Physical constants
    static constexpr double kB = 8.617333262145e-5; // eV K^-1 (Boltzmann constant)
    static constexpr double u = 931.49410242e6; // eV c^-2 (Atomic mass unit)

};