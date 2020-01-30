#include <algorithm>
#include <gsl/gsl_spline.h>
#include <memory>

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

	double result = gsl_spline_eval_integ(inter, x[0], x[x.size()-1], acc);

	gsl_spline_free(inter);
	gsl_interp_accel_free(acc);

	return result;
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

double Integrator::riemann_2d(const vector<double> &x, const vector<double> &y, const vector<vector<double>> &z) const {
	
	double integral_corners = 	  (x[1]-x[0])*(y[1]-y[0])*z[0][0]
								+ (x[1]-x[0])*(y[y.size()-1] - y[y.size()-2])*z[0][y.size()-1]
								+ (x[x.size()-1]-x[x.size()-2])*(y[1]-y[0])*z[x.size()-1][0]
								+ (x[x.size()-1]-x[x.size()-2])*(y[y.size()-1] - y[y.size()-2])*z[x.size()-1][y.size()-1];

	double integral_sides = 0.;
	for(size_t i = 1; i < x.size() - 1; ++i){
		integral_sides += (z[i][0]*(y[1]-y[0]) + z[i][y.size()-1]*(y[y.size()-1]-y[y.size()-2]))*(x[i]-x[i-1]);
	}
	for(size_t i = 1; i < y.size() - 1; ++i){
		integral_sides += (z[0][i]*(x[1]-x[0]) + z[x.size()-1][i]*(x[x.size()-1]-x[x.size()-2]))*(y[i]-y[i-1]);
	}

	double integral_area = 0.;
	for(size_t i = 1; i < x.size() - 1; ++i){
		for(size_t j = 1; j < y.size() - 1; ++j){
			integral_area += z[i][j]*(x[i]-x[i-1])*(y[j]-y[j-1]);
		}
	}

	return 0.25*integral_corners + 0.5*integral_sides + integral_area;
}