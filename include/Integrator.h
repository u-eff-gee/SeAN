/*    
    This file is part of SeAN.

    SeAN is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SeAN is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SeAN.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef INTEGRATOR_H
#define INTEGRATOR_H 1

#include <utility>
#include <vector>

using std::pair;
using std::vector;

class Integrator{

public:
	Integrator(){};
	~Integrator(){};

	// Trapezoidal rule for 1D and 2D histograms with equidistant bins
	double trapezoidal_rule(const vector<double> &bins, const vector<double> &histogram);
	double trapezoidal_rule2D(const vector<double> &bins1, const vector<double> &bins2, const vector<vector<double> > &histogram);
	// Simpson's rule (higher-order trapezoidal rule) for 2D histograms with equidistant bins
	double simpsons_rule2D(const vector<double> &bins1, const vector<double> &bins2, const vector<vector<double> > &histogram);

	// Darboux/Riemann integral for 1D and 2D functions 
	// Calculates the lower and upper sum of a function.
	// In contrast to the general definition of the lower and upper Darboux sum,
	// which is the infimum/supremum of a finite continuous interval, this implementation simply
	// takes the minimum/maximum of the boundaries of the interval.
	pair<double, double> darboux(const vector<double> &bins, const vector<double> &histogram);
	pair<double, double> darboux2D(const vector<double> &bins1, const vector<double> &bins2, const vector<vector<double> > &histogram);
};

#endif
