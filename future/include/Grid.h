#pragma once

#include<vector>

using std::vector;

class Grid{

public:
    Grid() = default;
    ~Grid() = default;

    vector<double> breit_wigner(const double x_0, const double gamma,
        const double x_min, const double x_max, const size_t n) const;
    double breit_wigner_coverage(const double x_0, const double gamma,
        const double x_min, const double x_max) const;

    vector<double> normal(const double mu, const double sigma,
        const double x_min, const double x_max, const size_t n) const;
    double normal_coverage(const double mu, const double sigma,
        const double x_min, const double x_max) const;

    void strictly_increasing(vector<double> &x) const;
};