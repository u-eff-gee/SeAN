#ifndef INTEGRATOR_H
#define INTEGRATOR_H 1

#include <vector>

using std::vector;

class Integrator{

public:
	Integrator(){};
	~Integrator(){};

	// Simple Riemann integrator for a 2D histogram with equidistant bins
	double integrate2DHistogram(const vector<double> &bins1, const vector<double> &bins2, const vector<vector<double> > &histogram);
};

#endif
