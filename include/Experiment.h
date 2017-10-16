#ifndef	EXPERIMENT_H 
#define EXPERIMENT_H 1

#include <vector>

#include "Target.h"

using std::vector;

class Experiment{

private:
	vector<double> energy_bins;
	vector<Target*> targets;
	vector<double> beamParams;

	string ifname;
	string beam_ID;

	double emin;
	double emax;

	unsigned int nBins_e;
	unsigned int nBins_z;

public:
	Experiment(): nBins_e(0), nBins_z(0){};
	~Experiment(){};

	// Functions to manage the calculation process
	void readInputFile (const char* filename);
	void crossSections(bool plot, bool write, bool exact);
	void transmission(bool plot, bool write);
	void print();

	// Functions to return private members
	Target* getTarget(unsigned int i){ return targets[i]; };

private:
	void createEnergyBins(double emin, double emax);
};

#endif
