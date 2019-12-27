#include <vector>

#include "Math/QuantFuncMathCore.h"

#include "Constants.h"
#include "Integrator.h"
#include "TestUtilities.h"

using std::vector;

double normal_pdf_2d(const double x, const double y, const double sigma_x, const double sigma_y) const {
    return Constants::inverse_twopi;
}

int main(){

    Integrator inte;
    TestUtilities test;
    // Test case for 2D integration:
    // Integrate a 2D normal distribution on a non-equidistant grid.
    const size_t n = 500;
    const double sigma{1.}; 

    vector<double> x(n), y(n);
    vector<vector<double>> z(n, vector<double>(n));

    for(size_t i = 0; i < n; ++i){
        for(size_t j = 0; j < n; ++j){
            x[i] = ROOT::Math::normal_quantile(1., sigma);
        }
    }
}