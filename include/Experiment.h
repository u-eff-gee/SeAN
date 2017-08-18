#ifndef	EXPERIMENT_H 
#define EXPERIMENT_H 1

#include <vector>

#include "Target.h"

using std::vector;

class Experiment{

private:
	double energy_bins[NBINS] = {0.};

	vector<Target*> targets;
	vector<double> beamParams;

	string ifname;
	string beam_ID;

	double emin;
	double emax;

public:
	Experiment(){};
	~Experiment(){};

	// Functions to manage the calculation process
	void readInputFile(const char* filename);
	void testIntegration(bool plot);
	void crossSections(bool plot, bool output);
	void transmission(bool plot, bool output);
	void print();
	void runTests();

	// Functions to return private members
	Target* getTarget(unsigned int i){ return targets[i]; };

private:
	void createEnergyBins(double emin, double emax);
};

#endif
