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

	Settings settings;

public:
	Experiment(Settings &s);
	~Experiment(){};

	// Functions to manage the calculation process
	void initialize();
	void crossSections();
	void transmission();

private:
	void createEnergyBins(double emin, double emax);
	void createTargets();
};

#endif
