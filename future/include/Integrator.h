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
	/// The trapezoidal rule linearly interpolates between the sampled points \f$x_i\f$ and \f$y_i\f$:
	/// \f[
	/// 	\int_{x_i}^{x_{i+1}} f(x) \mathrm{d} x \approx (x_{i+1} - x_i) \frac{f(x_{i+1}) + f(x_i)}{2}
	/// \f]
	double trapezoidal_rule(const vector<double> &x, const vector<double> &y) const;
	
	/// \brief (Cubic) Spline integration
	/// \param x Vector of x values
	/// \param y Vector of y values
	///
	/// The spline interpolation interpolates the sampled points \f$x_i\f$ and \f$y_i\f$ by a cubic spline.
	/// \cite Steffen1990
	double spline(const vector<double> &x, const vector<double> &y) const;

	pair<double, double> darboux(const vector<double> &x, const vector<double> &y) const;

};
