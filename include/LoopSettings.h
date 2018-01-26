#ifndef LOOPSETTINGS_H
#define LOOPSETTINGS_H 1

#include <string>
#include <vector>

#include "Settings.h"

using std::string;
using std::vector;

// Compactly process and store all necessary information to loop over SeAN input. Compared to Settings.h, all vectors have one additional dimension to store multiple values for a single parameter.
struct LoopSettings{

	// Settings for Experiment
	vector<double> emin;
	vector<double> emax;
	vector<double> nbins_e;
	vector<double> nbins_z;

	incidentBeamModel incidentBeam;
	vector< vector<double> > incidentBeamParams;
	string incidentBeamFile;

	// Settings for Targets
	vector<string> targetNames;
	vector<vector<vector<double> > > energy;
	vector<vector<vector<double> > > gamma0;
	vector<vector<vector<double> > > gamma;
	vector<vector<double> > ji;
	vector<vector<vector<double> > > jj;

	vector<vDistModel> vDist;
	vector< vector< vector<double> > > vDistParams;
	vector<string> vDistFile;

	vector<vector<double> > mass;

	vector<mAttModel> mAtt;
	vector<vector<vector<double> > > mAttParams;
	vector<vector<string> > mAttFile;

	vector<vector<double> > thickness;
	vector<vector<double> > velocity;
};

#endif
