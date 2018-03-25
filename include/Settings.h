#ifndef SETTINGS_H
#define SETTINGS_H 1

#include <string>
#include <vector>

using std::string;
using std::vector;

// Enum for different doppler-broadening models
// arb_vdist: arbitrary, point-wise defined velocity distribution
// arb_cs: arbitrary, point-wise defined velocity distribution
// zero: T = 0 -> no Doppler broadening of the cross section
// mb: Maxwell-Boltzmann distribution with effective temperature
// mba: Maxwell-Boltzmann distribution with effective temperature + use the approximation that Gamma is much smaller than the Doppler width
// phdos: Calculate the cross section from the eigenmodes (phonons) of the material
enum class dopplerModel{arb_vdist, arb_cs, zero, mb, mba, mbd, mbad, phdos};

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
        string inputfile = "";
	string outputfile = "";
	int verbosity = 2;
	bool recoil = false;
	bool output = false;
	bool exact = false;
	bool plot = false;
	bool write = false;

	// Settings for Experiment
	double emin = 0.;
	double emax = 0.;
	unsigned int nbins_e = 0;
	unsigned int nbins_z = 0;

	incidentBeamModel incidentBeam = incidentBeamModel::constant;
	vector<double> incidentBeamParams;
	string incidentBeamFile = "";

	// Settings for Targets
	vector<string> targetNames;
	vector<vector<double> > energy;
	vector<vector<double> > gamma0;
	vector<vector<double> > gamma;
	vector<double> ji;
	vector<vector<double> > jj;

	vector<dopplerModel> dopplerBroadening;
	vector< vector<double> > dopplerParams;
	vector<string> energyBinFile;
	vector<string> crosssectionFile;
	vector<string> velocityBinFile;
	vector<string> vDistFile;
	vector<string> omegaFile;
	vector<string> polarizationFile;
	vector<string> momentumFile;

	vector<double> mass;

	vector<mAttModel> mAtt;
	vector<vector<double> > mAttParams;
	vector<string> mAttFile;

	vector<double> thickness;
	vector<double> velocity;

	// Methods to print settings
	void print();
	void printOptions();
	void printExperiment();
	void printTarget(unsigned int i);

	// Method to write settings to file
	void write_output(unsigned int n_setting) const;
	void writeOptions(unsigned int n_setting) const;
	void writeExperiment() const;
	void writeTarget(unsigned int i) const;
};

#endif
