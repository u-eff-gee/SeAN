/*    
    This file is part of SeAN.

    SeAN is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SeAN is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SeAN.  If not, see <http://www.gnu.org/licenses/>.
*/


#pragma once 

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
	string outputfile = "out";
	int verbosity = 2;
	bool direct = false;
	bool recoil = false;
	bool output = false;
	bool exact = false;
	bool plot = false;
	bool status = false;
	bool write = false;
	bool write_all = false;
	bool uncertainty = false;
	bool multi = false;

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
	//vector<string> polarizationFile;
	//vector<string> momentumFile;

	vector<double> mass;

	vector<mAttModel> mAtt;
	vector<vector<double> > mAttParams;
	vector<string> mAttFile;

	vector<double> thickness;
	vector<double> velocity;

	// Methods to create output strings
	string bool_string(const bool b) const;
	string option_string() const;
	string experiment_string() const;
	string target_string(unsigned int i) const;

	// Methods to print settings
	void print();
	void printOptions();
	void printExperiment();
	void printTarget(unsigned int i);

	// Method to write settings to file
	void write_output() const;
	void writeOptions() const;
	void writeExperiment() const;
	void writeTarget(unsigned int i) const;
};

#endif
