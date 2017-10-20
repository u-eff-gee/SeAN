#ifndef	EXPERIMENT_H 
#define EXPERIMENT_H 1

#include <vector>

#include "Target.h"
#include "Settings.h"

using std::vector;

class Experiment{

private:
	vector<double> energy_bins;
	vector<Target*> targets;
	vector<double> beamParams;

	string ifname;
	string beam_ID;

	Settings settings;

public:
	Experiment(Settings &s);
	~Experiment(){};

	// Functions to manage the calculation process
	void readInputFile (const char* filename);
	void crossSections();
	void transmission();
	void print();

	// Functions to return private members
	Target* getTarget(unsigned int i){ return targets[i]; };

private:
	void createEnergyBins(double emin, double emax);
};

#endif
