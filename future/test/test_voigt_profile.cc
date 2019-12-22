#include <cmath>

#include "Constants.h"
#include "TestUtilities.h"
#include "VoigtProfile.h"

using std::sqrt;

int main(){
    TestUtilities test;

    VoigtProfile voigt;

    // Evaluate the Voigt profile at x = 0 for the two limits
    //  I) sigma = 1, gamma = 0, i.e. a normal distribution, which should evaluate to 1/sqrt(2*pi)
    test.test_closeness(voigt.pdf(0., 1., 0.),
                            Constants::inverse_sqrt_twopi,
                            test.num_tol_rel);
    // II) sigma = 0, gamma = 1, i.e. a Breit-Wigner distribution, which should evaluate to 2/pi
    test.test_closeness(voigt.pdf(0., 0., 1.),
                            2*Constants::inverse_pi,
                            test.num_tol_rel);
}