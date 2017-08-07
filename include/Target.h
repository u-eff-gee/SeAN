#ifndef TARGET_H 
#define TARGET_H 1

#include <string>
#include <vector>

#include "Config.h"
#include "CrossSection.h"
#include "Absorption.h"

using std::string;
using std::vector;

class Target{
private:
	vector< vector<double> > crosssection_bins;
	vector< vector<double> > velocity_bins;
	vector< vector<double> > vdist_bins;

	vector<double> e0_list;
	vector<double> gamma0_list;
	vector<double> gamma_list;
	vector<double> jj_list;
	vector<double> vDistParams;
	vector<double> vdist_norm;

	CrossSection *crossSection;
	Absorption *absorption;

	double incident_beam_bins[NBINS] = {1.};
	double dopplercs_bins[NBINS] = {0.};
	double massattenuation_bins[NBINS] = {0.};
	double transmitted_beam_bins[NBINS] = {0.};

	string target_name;
	string vDist_ID;
	string massAttenuation_ID;
	double j0;
	double mass;
	double z;

	int target_number;

public:	
	Target(string name, int number){
		target_name = name;
		target_number = number;

		crossSection = new CrossSection();
		absorption = new Absorption();
	};
	
	~Target(){
		delete crossSection;
	};

	void addEnergy(double e){ e0_list.push_back(e); };
	void addGamma0(double g0){ gamma0_list.push_back(g0); };
	void addGamma(double g){ gamma_list.push_back(g); };
	void addJJ(double j){ jj_list.push_back(j); };
	void setJ0(double j){ j0 = j; };
	void setMass(double m){ mass = m; };
	void setZ(double zz){ z = zz; };

	void setVDist(string vDist){vDist_ID = vDist; };
	void addVDistParameter(double p){ vDistParams.push_back(p); };

	void setMassAttenuation(string mAtt){ massAttenuation_ID = mAtt; };

	void print();

	// Functions to calculate histograms
	void calculateCrossSection(double (&energy_bins)[NBINS]);
	void calculateVelocityDistribution(double (&energy_bins)[NBINS]);
	void calculateDopplerShift(double (&energy_bins)[NBINS]);
	void calculateMassAttenuation(double (&energy_bins)[NBINS]);

	// Functions to plot histograms
	void plotCrossSection(double (&energy_bins)[NBINS]);
	void plotVelocityDistribution();
	void plotDopplerShift(double (&energy_bins)[NBINS]);
	void plotMassAttenuation(double (&energy_bins)[NBINS]);
	void plotMu();

private:
	double normalizeVDist(int i);
};

#endif
