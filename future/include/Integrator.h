#pragma once

#include <utility>
#include <vector>

using std::pair;
using std::vector;

class Integrator{

public:

	Integrator(){};
	~Integrator(){};

	double trapezoidal_rule(const vector<double> &x, const vector<double> &y) const;

	pair<double, double> darboux(const vector<double> &x, const vector<double> &y) const;

};