#ifndef	EXPERIMENT_H 
#define EXPERIMENT_H 1

#include <vector>

#include "Target.h"

using std::vector;

class Experiment{

private:
	double energy_bins[NBINS_E] = {0.};
	double velocity_bins[NBINS_V] = {0.};

	vector<Target*> targets;
	vector<double> beamParams;

	string ifname;
	string beam_ID;

	double emin;
	double emax;

public:
	Experiment(){};
	~Experiment(){};

	void readInputFile(const char* filename);
	void crossSections();
	void print();

private:
	void createEnergyBins(double emin, double emax);
};

#endif
