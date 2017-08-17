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

	vector<double> e0_at_rest_list;
	vector<double> e0_list;
	vector<double> gamma0_list;
	vector<double> gamma_list;
	vector<double> jj_list;
	vector<double> vDistParams;
	vector<double> vdist_norm;

	CrossSection *crossSection;
	Absorption *absorption;

	double photon_flux_density_bins[NBINS][NBINS_Z] = {{0.}, {0.}};
	double resonance_absorption_density_bins[NBINS][NBINS_Z] = {{0.}, {0.}};

	double incident_beam_bins[NBINS] = {1.};
	double dopplercs_bins[NBINS] = {0.};
	double massattenuation_bins[NBINS] = {0.};
	double z_bins[NBINS_Z] = {0.};
	double transmitted_beam_bins[NBINS] = {0.};


	string target_name;
	string vDist_ID;
	string massAttenuation_ID;
	double j0;
	double mass;
	double z;
	double vz;

	unsigned int target_number;

public:	
	Target(string name, unsigned int number){
		target_name = name;
		target_number = number;

		crossSection = new CrossSection();
		absorption = new Absorption();
	};
	
	~Target(){
		delete crossSection;
	};

	// Functions to set private member variables
	void addEnergy(double e){ e0_list.push_back(e); };
	void addEnergyAtRest(double e){ e0_at_rest_list.push_back(e); };
	void addGamma0(double g0){ gamma0_list.push_back(g0); };
	void addGamma(double g){ gamma_list.push_back(g); };
	void addJJ(double j){ jj_list.push_back(j); };
	void setJ0(double j){ j0 = j; };
	void setMass(double m){ mass = m; };
	void setZ(double zz){ z = zz; };
	void setVZ(double vvz){ vz = vvz; };

	void setVDist(string vDist){vDist_ID = vDist; };
	void addVDistParameter(double p){ vDistParams.push_back(p); };

	void setMassAttenuation(string mAtt){ massAttenuation_ID = mAtt; };

	// Function to modify resonance energies due to Doppler shift
	void boost();

	// Functions to return private member variables
	string getName(){ return target_name; };

	// Function to print information to the command line
	void print();

	// Functions to calculate histograms
	void calculateCrossSection(double (&energy_bins)[NBINS]);
	void calculateVelocityDistribution(double (&energy_bins)[NBINS]);
	void calculateDopplerShift(double (&energy_bins)[NBINS]);
	void calculateMassAttenuation(double (&energy_bins)[NBINS]);
	void calculateIncidentBeam(double (&energy_bins)[NBINS], string beam_ID, vector<double> beamParams);
	void calculateTransmittedBeam();
	void calculateZBins();
	void calculateZBins(double z0, double z1);
	void calculatePhotonFluxDensity();
	void calculateResonanceAbsorptionDensity();
	void calculateAbsorption(double (&energy_bins)[NBINS]);

	// Functions to plot histograms
	void plotCrossSection(double (&energy_bins)[NBINS]);
	void plotVelocityDistribution();
	void plotDopplerShift(double (&energy_bins)[NBINS]);
	void plotMassAttenuation(double (&energy_bins)[NBINS]);
	void plotMu();
	void plotPhotonFluxDensity(double (&energy_bins)[NBINS]);
	void plotTestIntegration(double (&energy_bins)[NBINS]);
	void plotResonanceAbsorptionDensity(double (&energy_bins)[NBINS]);

	// Function to set the incident beam and return the transmitted beam
	double& getTransmittedBeam(){ return transmitted_beam_bins[0]; };
	void setIncidentBeam(double &trans_beam_bins);

	// Functions to integrate over 2D histograms
	double integrateEZHistogram(double (&energy_bins)[NBINS], double (&z_bins)[NBINS_Z], double (&ezhist)[NBINS][NBINS_Z]);
	double integrateEEHistogram(double (&energy_bins)[NBINS], double (&eehist)[NBINS][NBINS]);

	// Functions for testing
	void testIntegration(double emin, double emax, double (&energy_bins)[NBINS], vector<double> beamParams);

private:
	double normalizeVDist(unsigned int i);
};

#endif
