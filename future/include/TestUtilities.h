#pragma once

#include <cmath>
#include <iostream>

using std::cout;
using std::endl;
using std::fabs;

class TestUtilities{

public:
    TestUtilities() = default;
    ~TestUtilities() = default;

    template<typename T, typename U>
    inline void test_equality(const T t, const U u) const {
        if(t != u){
            cout << __FILE__ << ": " << __FUNCTION__ << "( " << t << ", " << u << " )" << endl;
            throw;
        }
    }

    inline void test_closeness(const double x, const double y, const double epsilon_relative) const {
        if(fabs((x-y)/(0.5*(x+y))) > epsilon_relative){
            cout << __FILE__ << ": " << __FUNCTION__ << "( " << x << ", " << y << ", " << epsilon_relative << " )" << endl;
            throw;
        }
    }

    inline void test_larger_than(const double x, const double y) const {
        if(y > x){
            cout << __FILE__ << ": " << __FUNCTION__ << "( " << x << ", " << y << " )" << endl;
            throw;            
        }
    }

    static constexpr double num_tol_rel = 1e-5; // Standard required numerical precision
    static constexpr double num_tol_rel_ene = 1e-7; // Numerical precision for energies.
        // Due to the factor of 10^6 between eV and MeV, this value is smaller than 
        // num_tol_rel to be able to perform tests on the sub-eV level
    static constexpr double num_tol_rel_int = 1e-2; // Numerical precision for integrations.
        // The integrations are strongly dependent on the chosen range and the number 
        // of sampled points. The low threshold above is used for integrations of 
        // Breit-Wigner- and normal distributions with the same settings. For this reason,
        // it is considerably larger than num_tol_rel
};
