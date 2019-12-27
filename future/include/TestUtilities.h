#pragma once

#include <cmath>
#include <iostream>

using std::cout;
using std::endl;
using std::fabs;

/// \brief Utility functions and constants used for testing.
///
/// All methods with the naming pattern `test_*` will check whether a certain condition is fulfilled.
/// If this is not the case, an unspecified error will be thrown by involing `throw`.
/// This will cause `ctest` to report a failed test.
class TestUtilities{

public:
    /// \brief Default constructor.
    TestUtilities() = default;
    /// \brief Default destructor.
    ~TestUtilities() = default;

    /// \brief Function to test the equality of two values.
    /// \param t First value.
    /// \param u Second value.
    ///
    /// A comparison via `=` between the two types `T` and `U` must be possible in order for this function to work.
    template<typename T, typename U>
    inline void test_equality(const T t, const U u) const {
        if(t != u){
            cout << __FILE__ << ": " << __FUNCTION__ << "( " << t << ", " << u << " )" << endl;
            throw;
        }
    }

    /// \brief Function to test whether two numbers are close to each other.
    /// \param x First number.
    /// \param y Second number.
    /// \param epsilon_relative Maximum allowed relative deviation.
    ///
    /// This function tests whether
    /// \f[
    ///     \left| \frac{x-y}{\left(x+y\right)/2} \right| < \epsilon,
    /// \f]
    /// where \f$\epsilon\f$ is the maximum allowed relative deviation.
    inline void test_closeness(const double x, const double y, const double epsilon_relative) const {
        if(fabs((x-y)/(0.5*(x+y))) > epsilon_relative){
            cout << __FILE__ << ": " << __FUNCTION__ << "( " << x << ", " << y << ", " << epsilon_relative << " )" << endl;
            throw;
        }
    }

    /// \brief Function to test whether the first number is larger than the second one.
    /// \param x First number.
    /// \param y Second number.
    inline void test_larger_than(const double x, const double y) const {
        if(y > x){
            cout << __FILE__ << ": " << __FUNCTION__ << "( " << x << ", " << y << " )" << endl;
            throw;            
        }
    }

    static constexpr double num_tol_rel = 1e-5; ///< Standard required numerical precision for tests.
    /// \brief Numerical precision for energies.
    /// 
    /// Due to the factor of 10^6 between eV and MeV, this value is smaller than 
    /// `Test::num_tol_rel` to be able to perform tests on the sub-eV level.
    static constexpr double num_tol_rel_ene = 1e-7;
    /// \brief Numerical precision for integrations.
    ///
    /// The integrations are strongly dependent on the chosen range and the number 
    /// of sampled points. The low threshold below is used for integrations of 
    /// Breit-Wigner- and normal distributions with the same settings. For this reason,
    /// it is considerably larger than `Test::num_tol_rel`.
    static constexpr double num_tol_rel_int = 1e-2;
};
