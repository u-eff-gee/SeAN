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


#include "Integrator.h"

double Integrator::trapezoidal_rule(const vector<double> &bins, const vector<double> &histogram){
	// Bin size
	double bin_size = bins[1] - bins[0];

	double integral = 0.;

	#pragma omp parallel for reduction (+:integral)
	for(unsigned int i = 0; i < bins.size() - 1; ++i){
		integral += (histogram[i] + histogram[i+1]);
	}

	return 0.5*bin_size*integral;
}

double Integrator::integrate2DHistogram(const vector<double> &bins1, const vector<double> &bins2, const vector<vector<double> > &histogram){

	// Area of a bin in 2D plane
	double bin_area = (bins1[1] - bins1[0])*(bins2[1] - bins2[0]);

	// Crude implementation of the Riemann integral as a sum of bin contents times their dimension in energy- and z-direction. Since there are only settings.nbins_e-1 spaces between settings.nbins_e bins, leave out the last bin in each loop.
	double integral = 0.;

	#pragma omp parallel for reduction (+:integral)
	for(unsigned int i = 0; i < bins1.size() - 1; ++i){
		for(unsigned int j = 0; j < bins2.size() - 1; ++j){
			integral += histogram[i][j]; 
		}
	}

	return bin_area*integral;
}
