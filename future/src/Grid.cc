#include "Math/ProbFuncMathCore.h"
#include "Math/QuantFuncMathCore.h"

#include "Grid.h"

vector<double> Grid::breit_wigner(const double x_0, const double gamma,
    const double x_min, const double x_max, const size_t n_x) const {

    vector<double> grid(n_x, 0.);
    grid[0] = x_min;
    grid[n_x-1] = x_max;

    double quantile_min = ROOT::Math::breitwigner_cdf(x_min, gamma, x_0);
    double quantile_max = ROOT::Math::breitwigner_cdf(x_max, gamma, x_0);
    double quantile_inc = (quantile_max-quantile_min)/(n_x-1);

    for(size_t i = 1; i < n_x-1; ++i){
        grid[i] = ROOT::Math::breitwigner_quantile(quantile_min+i*quantile_inc, gamma) + x_0;
    }

    return grid;
}

double Grid::breit_wigner_coverage(const double x_0, const double gamma,
    const double x_min, const double x_max) const {
        return ROOT::Math::breitwigner_cdf(x_max, gamma, x_0)
            -ROOT::Math::breitwigner_cdf(x_min, gamma, x_0);
}

vector<double> Grid::normal(const double mu, const double sigma,
    const double x_min, const double x_max, const size_t n_x) const {

    vector<double> grid(n_x, 0.);
    grid[0] = x_min;
    grid[n_x-1] = x_max;

    double quantile_min = ROOT::Math::normal_cdf(x_min, sigma, mu);
    double quantile_max = ROOT::Math::normal_cdf(x_max, sigma, mu);
    double quantile_inc = (quantile_max-quantile_min)/(n_x-1);

    for(size_t i = 1; i < n_x-1; ++i){
        grid[i] = ROOT::Math::normal_quantile(quantile_min+i*quantile_inc, sigma) + mu;
    }

    return grid;
}

double Grid::normal_coverage(const double mu, const double sigma,
    const double x_min, const double x_max) const {
        return ROOT::Math::normal_cdf(x_max, sigma, mu)
            -ROOT::Math::normal_cdf(x_min, sigma, mu);
}