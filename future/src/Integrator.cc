#include <algorithm>
#include <gsl/gsl_spline.h>

#include "Integrator.h"

double Integrator::trapezoidal_rule(const vector<double> &x, const vector<double> &y) const {

	double integral = 0.;

	for(unsigned int i = 0; i < x.size() - 1; ++i){
		integral += (y[i+1] + y[i])*(x[i+1]-x[i]);
	}

	return 0.5*integral;
}
	
double Integrator::spline(const vector<double> &x, const vector<double> &y) const {
	
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	gsl_spline *inter = gsl_spline_alloc(gsl_interp_steffen, x.size());
	gsl_spline_init(inter, &x[0], &y[0], x.size());

	return gsl_spline_eval_integ(inter, x[0], x[x.size()-1], acc);
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
