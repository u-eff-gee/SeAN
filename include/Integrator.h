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

#include <vector>

using std::vector;

class Integrator{

public:
	Integrator(){};
	~Integrator(){};

	// Trapezoidal rule for 1D histograms with equidistant bins
	double trapezoidal_rule(const vector<double> &bins, const vector<double> &histogram);
	// Simple Riemann integrator for a 2D histogram with equidistant bins
	double integrate2DHistogram(const vector<double> &bins1, const vector<double> &bins2, const vector<vector<double> > &histogram);
};

#endif
