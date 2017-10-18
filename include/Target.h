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
	vector<double> z_bins;

	vector<double> e0_at_rest_list;
	vector<double> e0_list;
	vector<double> gamma0_list;
	vector<double> gamma_list;
	vector<double> jj_list;
	vector<double> vDistParams;
	vector<double> vdist_norm;
	vector<unsigned int> vdist_centroid;

	CrossSection *crossSection;
	Absorption *absorption;

	vector<vector <double> > photon_flux_density_bins;
	vector<vector <double> > resonance_absorption_density_bins;

	vector<double> incident_beam_bins;
	vector<double> dopplercs_bins;
	vector<double> massattenuation_bins;
	vector<double> transmitted_beam_bins;


	string target_name;
	string vDist_ID;
	string massAttenuation_ID;
	double j0;
	double mass;
	double z;
	double vz;

	unsigned int target_number;

public:	
	Target(unsigned int ne, unsigned int nz, string name, unsigned int number) : z_bins(nz, 0.), incident_beam_bins(ne, 0.), dopplercs_bins(ne, 0.), massattenuation_bins(ne, 0.), transmitted_beam_bins(ne, 0.){
		target_name = name;
		target_number = number;

		crossSection = new CrossSection();
		absorption = new Absorption();

		photon_flux_density_bins.reserve(ne*nz);
		resonance_absorption_density_bins.reserve(ne*nz);

		for(unsigned int i = 0; i < nz; ++i){
			photon_flux_density_bins.push_back(vector<double>(ne, 0.));
			resonance_absorption_density_bins.push_back(vector<double>(ne, 0.));
		}
	};
	
	~Target(){
		delete crossSection;
		delete absorption;

		delete &crosssection_bins;
		delete &velocity_bins;
		delete &vdist_bins;
		delete &z_bins;
		delete &e0_at_rest_list;
		delete &e0_list;
		delete &gamma0_list;
		delete &gamma_list;
		delete &jj_list;
		delete &vDistParams;
		delete &vdist_norm;
		delete &vdist_centroid;

		delete &photon_flux_density_bins;
		delete &resonance_absorption_density_bins;

		delete &incident_beam_bins;
		delete &dopplercs_bins;
		delete &massattenuation_bins;
		delete &transmitted_beam_bins;
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

	// Function to read the nuclear mass from the AME table
	void readAME(string isotope);

	// Functions to return private member variables
	string getName(){ return target_name; };

	// Function to print information to the command line
	void print();

	// Functions to calculate histograms
	void calculateCrossSection(vector<double> &energy_bins);
	void calculateVelocityDistribution(vector<double> &energy_bins);
	void calculateDopplerShift(vector<double> &energy_bins);
	void calculateDopplerShiftFFT(vector<double> &energy_bins);
	void calculateMassAttenuation(vector<double> &energy_bins);
	void calculateIncidentBeam(vector<double> &energy_bins, string beam_ID, vector<double> beamParams);
	void calculateTransmittedBeam();
	void calculateZBins();
	void calculateZBins(double z0, double z1);
	void calculatePhotonFluxDensity();
	void calculateResonanceAbsorptionDensity();
	void calculateAbsorption(vector<double> &energy_bins);

	// Functions to plot histograms
	void plotCrossSection(vector<double> &energy_bins);
	void plotVelocityDistribution();
	void plotDopplerShift(vector<double> &energy_bins);
	void plotMassAttenuation(vector<double> &energy_bins);
	void plotMu();
	void plotPhotonFluxDensity(vector<double> &energy_bins);
	void plotTestIntegration(vector<double> &energy_bins);
	void plotResonanceAbsorptionDensity(vector<double> &energy_bins);

	// Functions to write histograms to file
	void write(vector<double> &energy_bins);
	void print1DVector(vector<double> &vec, string column, string filename);
	void print2DVector(vector<vector<double> > &vec, string column, string filename);

	// Function to set the incident beam and return the transmitted beam
	double& getTransmittedBeam(){ return transmitted_beam_bins[0]; };
	void setIncidentBeam(double &trans_beam_bins);

	// Functions to integrate over 2D histograms
	double integrateEZHistogram(vector<double> &energy_bins, vector<double> &z_bins, vector<vector<double> > &ezhist);
	double integrateEEHistogram(vector<double> &energy_bins, vector<vector<double> > &eehist);

	void vDistInfo(vector< vector<double> > &velocity_bins, vector< vector<double> > &vdist_bins, vector<double> &vdist_norm, vector<unsigned int> &vdist_centroid);
};

#endif
