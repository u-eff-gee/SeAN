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

};
