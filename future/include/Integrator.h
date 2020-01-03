#pragma once

#include <utility>
#include <vector>

using std::pair;
using std::vector;

/// \brief Class for the integration of functions
class Integrator{

public:

	/// \brief Default constructor
	Integrator() = default;
	/// \brief Default destructor
	~Integrator() = default;

	/// \brief Trapezoidal rule integration
	/// \param x Vector of x values
	/// \param y Vector of y values
	///
	/// The trapezoidal rule \cite TrapezoidalRule2019 linearly interpolates between the sampled points \f$x_i\f$ and \f$y_i\f$:
	/// \f[
	/// 	\int_{x_i}^{x_{i+1}} f(x) \mathrm{d} x \approx (x_{i+1} - x_i) \frac{f(x_{i+1}) + f(x_i)}{2}
	/// \f]
	double trapezoidal_rule(const vector<double> &x, const vector<double> &y) const;
	
	/// \brief (Cubic) Spline integration
	/// \param x Vector of x values
	/// \param y Vector of y values
	///
	/// The spline interpolation interpolates the sampled points \f$x_i\f$ and \f$y_i\f$ by a cubic spline which conserves their monotonicity \cite Steffen1990.
	/// After that, the cubic spline is integrated.
	/// This function uses the so-called 'Steffen interpolation' from the GNU Scientific Library (GSL) \cite Galassi2019.
	double spline(const vector<double> &x, const vector<double> &y) const;

	/// \brief Darboux integration (lower and upper sum)
	/// \param x Vector of x values
	/// \param y Vector of y values
	/// \returns Value of the integral as a `std::pair<double, double>`, where `pair.first` is the lower- and `pair.second` the upper sum.
	///
	/// The Darboux integral \cite DarbouxIntegral2019 uses a stepwise interpolation of the sampled points \f$x_i\f$ and \f$y_i\f$ to compute 
	/// a lower sum \f$L\f$ and an upper sum \f$U\f$:
	/// \f[
	/// 	\int_{x_i}^{x_{i+1}} f(x) \mathrm{d} x \approx L = (x_{i+1} - x_i) \mathrm{min} \left[ \left\{ f(x_i), f(x_{i+1}) \right\} \right]
	/// \f]
	/// \f[
	/// 	\int_{x_i}^{x_{i+1}} f(x) \mathrm{d} x \approx U = (x_{i+1} - x_i) \mathrm{max} \left[ \left\{ f(x_i), f(x_{i+1}) \right\} \right],
	/// \f]
	/// whose difference can be made arbitrarily small for a Darboux-integrable function.
	/// This function computes both \f$L\f$ and \f$U\f$ and returns them as a `std::pair<double, double>`.
	/// For a sufficiently fine grid, they can be seen as strict lower and upper limits of the trapezoidal rule and the spline interpolation.
	pair<double, double> darboux(const vector<double> &x, const vector<double> &y) const;

	/// \brief Two-dimensional Riemann integration
	/// \param x Vector of x values
	/// \param y Vector of y values
	/// \param z Vector of z values
	///
	/// This 2D integration algorithm, which is inspired by the 1D Riemann integral \cite RiemannIntegral2020, 
	/// approximates the 2D function by rectangular cuboids:
	/// \f[
	///		\int_{x_i}^{x_{i+1}} \int_{y_j}^{y_{j+1}} f(x, y) \mathrm{d} x \mathrm{d} y \approx f(x_i, y_i) \left( x_{i+1} - x_i \right) \left( y_{i+1} - y_i \right).
	/// \f]
	/// The integration over the entire domain \f$x \in \left[ x_\mathrm{start}, x_\mathrm{stop} \right]\f$ and 
	/// \f$y \in \left[ y_\mathrm{start}, y_\mathrm{stop} \right]\f$ is performed by sampling \f$f\f$ at generally
	/// non-equidistant points \f$x_i\f$ and \f$y_i\f$. The number of points in \f$x\f$- and \f$y\f$ direction is
	/// given by \f$n_x\f$ and \f$n_y\f$. Three distinct sums over the corners (\f$C\f$), the sides (\f$S\f$), and
	/// the inner area (\f$I\f$) of the rectangular domain are calculated:
	/// \f[
	/// 	C = f(x_0, y_0) \left( x_1 - x_0 \right) \left( y_1 - y_0 \right) 
	///			+ f(x_0, y_{n_y - 1}) \left( x_1 - x_0 \right) \left( y_{n_y - 1} - y_{n_y - 2} \right)
	/// \f]
	/// \f[
	/// 	+ f(x_{n_x - 1}, y_0) \left( x_{n_x - 1} - x_{n_x - 2} \right) \left( y_1 - y_0 \right)
	/// 	+ f(x_{n_x - 1}, y_{n_y - 1}) \left( x_{n_x - 1} - x_{n_x - 2} \right) \left( y_{n_y - 1} - y_{n_y - 2} \right)
	/// \f]
	/// \f[
	/// 	S = \sum_{i=1}^{n_y - 2} \left[ f(x_0, y_i) \left( x_1 - x_0 \right) 
	///			+ f(x_{n_x - 1}, y_i) \left( x_{n_x - 1} - x_{n_x - 2} \right) \right]
	///			\left( y_i - y_{i-1} \right)
	/// \f]
	/// \f[
	/// 	+ \sum_{i=1}^{n_x - 2} \left[ f(x_i, y_0) \left( y_1 - y_0 \right) 
	///			+ f(x_i, y_{n_y - 1}) \left( y_{n_y - 1} - y_{n_y - 2} \right) \right]
	///			\left( x_i - x_{i-1} \right)
	/// \f]
	/// \f[
	/// 	I = \sum_{i = 1}^{n_x - 2} \sum_{j = 1}^{n_y - 2} f(x_i, y_j) \left( x_i - x_{i-1} \right) \left( y_i - y_{i-1} \right),
	/// \f]
	/// yielding the total value of the integral as:
	/// \f[
	///		\int_{x_\mathrm{start}}^{x_\mathrm{stop}} \int_{y_\mathrm{start}}^{y_\mathrm{stop}} f(x, y) \mathrm{d} x \mathrm{d} y \approx \frac{1}{4} C + \frac{1}{2} S + I.
	/// \f]
	/// The fractions in front of \f$C\f$ and \f$S\f$ correct for the double counting of certain grid points.
	double riemann_2d(const vector<double> &x, const vector<double> &y, const vector<vector<double>> &z) const;
};
