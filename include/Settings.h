#ifndef SETTINGS_H
#define SETTINGS_H 1

#include <string>
#include <vector>

using std::string;
using std::vector;

// Enum for different velocity distribution models
// arb: arbitrary, point-wise defined distribution
// zero: T = 0 -> no Doppler broadening of the cross section
// mb: Maxwell-Boltzmann distribution with effective temperature
// mba: Maxwell-Boltzmann distribution with effective temperature + use the approximation that Gamma is much smaller than the Doppler width
enum class vDistModel{arb, zero, mb, mba};

// Enum for different mass attenuation models
// arb: arbitrary, point-wise defined mass attenuation
// nist: using tabulated data from NIST
enum class mAttModel{constant, arb, nist};

// Enum for incident beam intensity distribution
// constant: constant intensity over the whole energy range
// gauss: normal distribution
// arb: arbitrary, point-wise defined intensity distribution
enum class incidentBeamModel{constant, gauss, arb};

// Stores all the input information for SeAN
struct Settings{
	// Command-line options
        char *inputfile;
	bool exact = false;
	bool plot = false;
	bool write = false;
	bool sudowrite = false;

	// Settings for Experiment
	double emin;
	double emax;
	unsigned int nbins_e;
	unsigned int nbins_z;

	incidentBeamModel incidentBeam;
	vector<double> incidentBeamParams;
	string incidentBeamFile;

	// Settings for Targets
	vector<string> targetNames;
	vector<vector<double> > energy;
	vector<vector<double> > gamma0;
	vector<vector<double> > gamma;
	vector<vector<double> > ji;
	vector<vector<double> > jj;

	vector<vDistModel> vDist;
	vector<double> vDistParams;

	vector<double> mass;

	vector<mAttModel> mAtt;
	vector<double> mAttParams;
	string mAttFile;

	vector<double> thickness;
	vector<double> velocity;

	void printOptions();
	void printExperiment();
	void printTarget(unsigned int i);
};

#endif
