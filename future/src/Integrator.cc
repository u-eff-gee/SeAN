#include <algorithm>

#include "Integrator.h"

double Integrator::trapezoidal_rule(const vector<double> &x, const vector<double> &y) const {

	double integral = 0.;

	for(unsigned int i = 0; i < x.size() - 1; ++i){
		integral += (y[i+1] + y[i])*(x[i+1]-x[i]);
	}

	return 0.5*integral;
}

pair<double, double> Integrator::darboux(const vector<double> &x, const vector<double> &y) const {

	double lower_sum = 0.;
	double upper_sum = 0.;

	for(unsigned int i = 0; i < x.size() - 1; ++i){

		lower_sum += std::min(y[i+1], y[i])*(x[i+1]-x[i]);
		upper_sum += std::max(y[i+1], y[i])*(x[i+1]-x[i]);
	
    }

	return pair<double, double>(lower_sum, upper_sum);
}