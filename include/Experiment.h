#ifndef	EXPERIMENT_H 
#define EXPERIMENT_H 1

#include <vector>
#include <iostream>

#include "Target.h"
#include "Settings.h"

using std::vector;
using std::endl;
using std::cout;

class Experiment{

private:
	vector<double> energy_bins;
	vector<Target> targets;

	Settings settings;

public:
	Experiment(){};
	Experiment(Settings &s);
	~Experiment(){
//		for(unsigned int i = 0; i < targets.size(); ++i){
//			delete targets[i];
//		}
	};

	// Functions to manage the calculation process
	void initialize();
	void crossSections();
	void transmission();
	void resonant_scattering();

	// Functions for output
	void plot(unsigned int n_setting) ;
	void print_results();
	void write(unsigned int n_setting) ;
	void write_results(string outputfile, unsigned int n_setting) const;

private:
	void createEnergyBins(double emin, double emax);
	void createTargets();
};

#endif
