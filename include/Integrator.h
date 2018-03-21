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
