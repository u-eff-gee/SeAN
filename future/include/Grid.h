#pragma once

#include<vector>

using std::vector;

/// \brief Class for the generation of grids for the evaluation of probability distribution functions (PDFs).
///
/// For any PDF \f$f(x)\f$, a grid will be generated from the corresponding quantile function \f$F^{-1}\f$, which is the inverse of the cumulative distribution function (CDF):
/// \f[
///     F(x) = \int_{-\infty}^{x} f(x'),
/// \f]
/// \f[
///     F(x) \in [0, 1].
/// \f]
/// The grid generation is related to the inverse transform sampling method \cite InverseTransformSampling2019 for generating random numbers from an arbitrary distribution.
/// It is based on the fact that a set of points \f$x_i\f$, obtained from a number of \f$n\f$ uniformly distributed points \f$u_i\f$ by
/// \f[
///     x_i = F^{-1}(u_i)
/// \f]
/// will be distributed according to \f$f\f$.
/// Given a start point \f$x_\mathrm{min}\f$, an end point \f$x_\mathrm{max}\f$, and a number of required points \f$n\f$, a set of equidistant \f$u_i\f$ will be used:
/// \f[
///     u_0 = F(x_\mathrm{min}),
/// \f]
/// \f[
///     u_{n-1} = F(x_\mathrm{max}),
/// \f]
/// \f[
///     u_{i} = u_0 + i \frac{u_{n-1} - u_0}{n - 1}, ~~~~ 0 < i < n-1
/// \f]
/// The coverage of the PDF by the grid is given by
/// \f[
///     F(x_\mathrm{max}) - F(x_\mathrm{min})
/// \f]
class Grid{

public:
    /// \brief Default constructor.
    Grid() = default;
    /// \brief Default destructor.
    ~Grid() = default;

    /// \brief Generate a grid for a Breit-Wigner (Cauchy) distribution \cite CauchyDistribution2019.
    /// \param x_0 Centroid.
    /// \param gamma Scale parameter.
    /// \param x_min Start point of the grid.
    /// \param x_max End point of the grid.
    /// \param n Number of grid points.
    vector<double> breit_wigner(const double x_0, const double gamma,
        const double x_min, const double x_max, const size_t n) const;
    /// \brief Calculate the coverage of a Breit-Wigner distribution by the grid.
    /// \param x_0 Centroid.
    /// \param gamma Scale parameter.
    /// \param x_min Start point of the grid.
    /// \param x_max End point of the grid.
    double breit_wigner_coverage(const double x_0, const double gamma,
        const double x_min, const double x_max) const;
    
    /// \brief Generate a grid for a normal distribution \cite NormalDistribution2019.
    /// \param mu Mean value.
    /// \param sigma Standard deviation.
    /// \param x_min Start point of the grid.
    /// \param x_max End point of the grid.
    /// \param n Number of grid points.
    vector<double> normal(const double mu, const double sigma,
        const double x_min, const double x_max, const size_t n) const;
    
    /// \brief Calculate the coverage of a normal distribution by the grid
    /// \param mu Mean value.
    /// \param sigma Standard deviation.
    /// \param x_min Start point of the grid.
    /// \param x_max End point of the grid.
    /// \param n Number of grid points.
    double normal_coverage(const double mu, const double sigma,
        const double x_min, const double x_max) const;

    /// \brief Create a grid of strictly increasing values from a given set of values 
    /// \param x Set of x values which is potentially unsorted and may contain the same x value multiple times.
    ///
    /// First, the set of values \f$x_i\f$ is sorted to give a set of values \f$\tilde{x}_i\f$.
    /// After that, equal numbers are eliminated by interpolating linearly to the next higher different number, i.e.:
    /// \f[
    ///     \left\{ 1.0, 2.0, 2.0, 3.0 \right\} \to \left\{ 1.0, 2.0, 2.5, 3.0 \right\}
    /// \f]
    /// Interpolating to the next-higher number does not work if the set of equal numbers includes the last value in the list.
    /// In this special case, the interpolation of performed backwards, i.e.:
     /// \f[
    ///     \left\{ 1.0, 2.0, 2.0 \right\} \to \left\{ 1.0, 1.5, 2.0 \right\}
    /// \f]       
    void strictly_increasing(vector<double> &x) const;
};