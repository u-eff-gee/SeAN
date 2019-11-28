#include <vector>

#include "Integrator.h"
#include "TestUtilities.h"

using std::vector;

int main(){

    Integrator inte;
    TestUtilities test;

    // Simple test case:
    // Integration of the function y = x with non-equidistant values of x.
    vector<double> x{0., 0.05, 1.};
    double inte_trap = inte.trapezoidal_rule(x, x);
    test.is_close_relative(inte_trap, 0.5, test.num_tol_rel);

    pair<double, double> inte_darb = inte.darboux(x, x);
    test.is_equal(inte_darb.first, 0.05*0.95);
    test.is_equal(inte_darb.second, 0.05*0.05+0.95*1.);

}