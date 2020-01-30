#include <vector>

#include "Math/QuantFuncMathCore.h"

#include "Constants.h"
#include "Integrator.h"
#include "TestUtilities.h"

using std::vector;

double normal_pdf_2d(const double x, const double y, const double sigma_x, const double sigma_y) {
    const double inverse_sigma_x = 1./sigma_x;
    const double inverse_sigma_y = 1./sigma_y;

    return Constants::inverse_twopi*inverse_sigma_x*inverse_sigma_y
                *exp(-0.5*(x*x*inverse_sigma_x*inverse_sigma_x
                            +y*y*inverse_sigma_y*inverse_sigma_y));
}

int main(){

    Integrator inte;
    TestUtilities test;
    // Test case for 2D integration:
    // Integrate a 2D normal distribution on a non-equidistant grid.
    const size_t n = 1000;
    const double sigma_x{1.}, sigma_y{1.};
    const double quantile_min = 1e-5;
    const double quantile_inc = (1.-2*quantile_min)/(n-1);

    vector<double> x(n), y(n);
    vector<vector<double>> z(n, vector<double>(n));

    for(size_t i = 0; i < n; ++i){
            x[i] = ROOT::Math::normal_quantile(i*quantile_inc + quantile_min, sigma_x);
            y[i] = ROOT::Math::normal_quantile(i*quantile_inc + quantile_min, sigma_y);
    }

    for(size_t i = 0; i < n; ++i){
        for(size_t j = 0; j < n; ++j){
            z[i][j] = normal_pdf_2d(x[i], y[j], sigma_x, sigma_y);
        }
    }

    double inte_result = inte.riemann_2d(x, y, z);

    test.test_closeness(inte_result, 1., test.num_tol_rel_int);
}