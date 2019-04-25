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

#include <algorithm>
#include "Integrator.h"

using std::max;
using std::max_element;
using std::min;
using std::min_element;

double Integrator::trapezoidal_rule(const vector<double> &bins, const vector<double> &histogram){

	double bin_size = bins[1] - bins[0];

	double integral = 0.;

	#pragma omp parallel for reduction (+:integral)
	for(unsigned int i = 0; i < bins.size() - 1; ++i){
		integral += (histogram[i] + histogram[i+1]);
	}

	return 0.5*bin_size*integral;
}

double Integrator::trapezoidal_rule2D(const vector<double> &bins1, const vector<double> &bins2, const vector<vector<double> > &histogram){

	double bin_area = (bins1[1] - bins1[0])*(bins2[1] - bins2[0]);

	double integral_sides = 0.;
	double integral_inside = 0.;

	// Sum over corners
	double integral_corners = (histogram[0             ][0] + histogram[0             ][bins2.size()-1] +
				   histogram[bins1.size()-1][0] + histogram[bins1.size()-1][bins2.size()-1]);

	// Sum over y = 0 and y = (bins2.size()-1) axis
	#pragma omp parallel for reduction (+:integral_sides)
	for(unsigned int i = 1; i < bins1.size() - 1; ++i){
		integral_sides += histogram[i][0             ] +
				  histogram[i][bins2.size()-1];
	}
	// Sum over x = 0 and x = (bins1.size()-1) axis
	#pragma omp parallel for reduction (+:integral_sides)
	for(unsigned int i = 1; i < bins2.size() - 1; ++i){
		integral_sides += histogram[0             ][i] +
				  histogram[bins1.size()-1][i];
	}

	// Sum over inner area
	#pragma omp parallel for reduction (+:integral_inside)
	for(unsigned int i = 2; i < bins1.size() - 2; ++i){
		for(unsigned int j = 2; j < bins2.size() - 2; ++j){
			integral_inside += histogram[i][j];
		}
	}
	return bin_area*(0.25*integral_corners + 0.5*integral_sides + integral_inside);
}

double Integrator::simpsons_rule2D(const vector<double> &bins1, const vector<double> &bins2, const vector<vector<double> > &histogram){

	double bin_area = (bins1[1] - bins1[0])*(bins2[1] - bins2[0]);

	double integral_value_2_sides = 0.;
	double integral_value_4_sides = 0.;

	double integral_value_4_area = 0.;
	double integral_value_8_area = 0.;
	double integral_value_16_area = 0.;

	// Sum over corners
	double integral_corners = (histogram[0             ][0] + histogram[0             ][bins2.size()-1] +
				   histogram[bins1.size()-1][0] + histogram[bins1.size()-1][bins2.size()-1]);

	// Sum over sides
	// Sum over even bins on y = 0 and y = (bins2.size()-1) axis
	#pragma omp parallel for reduction (+:integral_value_2_sides)
	for(unsigned int i = 2; i < bins1.size() - 1; i += 2){
		integral_value_2_sides += histogram[i][0             ] +
				  		histogram[i][bins2.size()-1];
	}
	// Sum over odd bins on y = 0 and y = (bins2.size()-1) axis
	#pragma omp parallel for reduction (+:integral_value_4_sides)
	for(unsigned int i = 1; i < bins1.size() - 2; i += 2){
		integral_value_4_sides += histogram[i][0             ] +
				  		histogram[i][bins2.size()-1];
	}
	// Sum over even bins on x = 0 and x = (bins2.size()-1) axis
	#pragma omp parallel for reduction (+:integral_value_2_sides)
	for(unsigned int i = 2; i < bins2.size() - 1; i += 2){
		integral_value_2_sides += histogram[0             ][i] +
				  		histogram[bins1.size()-1][i];
	}
	// Sum over odd bins on x = 0 and x = (bins2.size()-1) axis
	#pragma omp parallel for reduction (+:integral_value_4_sides)
	for(unsigned int i = 1; i < bins2.size() - 2; i += 2){
		integral_value_4_sides += histogram[0             ][i] +
						histogram[bins1.size()-1][i];
	}

	// Sum over inner area
	// Sum over inner even-even bins
	#pragma omp parallel for reduction (+:integral_value_4_area)
	for(unsigned int i = 2; i < bins1.size()-1; i +=2){
		for(unsigned int j = 2; j < bins2.size()-1; j +=2){
			integral_value_4_area += histogram[i][j];
		}
	}
	// Sum over inner odd-odd bins
	#pragma omp parallel for reduction (+:integral_value_16_area)
	for(unsigned int i = 1; i < bins1.size()-2; i +=2){
		for(unsigned int j = 1; j < bins2.size()-2; j +=2){
			integral_value_16_area += histogram[i][j];
		}
	}
	// Sum over inner odd-even bins
	#pragma omp parallel for reduction (+:integral_value_8_area)
	for(unsigned int i = 1; i < bins1.size()-2; i +=2){
		for(unsigned int j = 2; j < bins2.size()-1; j +=2){
			integral_value_8_area += histogram[i][j] + histogram[i+1][j-1];
		}
	}

	return 1./9.*bin_area*(integral_corners + 
			       2.*integral_value_2_sides + 4.*integral_value_4_sides +
			       4.*integral_value_4_area  + 8.*integral_value_8_area +
			       16.*integral_value_16_area);
}

pair<double, double> Integrator::darboux(const vector<double> &bins, const vector<double> &histogram){

	double bin_size = bins[1] - bins[0];

	double lower_sum = 0.;
	double upper_sum = 0.;
	pair<double, double> mima;

	#pragma omp parallel for reduction(+:lower_sum) reduction(+:upper_sum)
	for(unsigned int i = 0; i < bins.size() - 1; ++i){
		// Taking the minimum and the maximum of the same pair of values seems
		// inefficient. In fact, there is the function std::minmax which evaluates
		// both minimum and maximum value at the same time, but this did not work
		// well with omp parallelism.
		// Sacrificed the usage of minmax for omp.
		lower_sum += min(histogram[i], histogram[i+1]);
		upper_sum += max(histogram[i], histogram[i+1]);
	}

	return pair<double, double>(bin_size*lower_sum, bin_size*upper_sum);
}

pair<double, double> Integrator::darboux2D(const vector<double> &bins1, const vector<double> &bins2, const vector<vector<double> > &histogram){

	double bin_area = (bins1[1] - bins1[0])*(bins2[1] - bins2[0]);

	double lower_sum = 0.;
	double upper_sum = 0.;

	#pragma omp parallel for reduction (+:lower_sum) reduction (+:upper_sum)
	for(unsigned int i = 0; i < bins1.size() - 1; ++i){
		for(unsigned int j = 0; j < bins2.size() - 1; ++j){
			lower_sum += min(min(histogram[i  ][j  ], histogram[i  ][j+1]),
					 min(histogram[i+1][j  ], histogram[i+1][j+1]));
			upper_sum += max(max(histogram[i  ][j  ], histogram[i  ][j+1]),
					 max(histogram[i+1][j  ], histogram[i+1][j+1]));
		}
	}

	return pair<double, double>(bin_area*lower_sum, bin_area*upper_sum);
}
