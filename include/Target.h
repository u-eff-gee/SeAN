#ifndef TARGET_H 
#define TARGET_H 1

#include <string>
#include <vector>

#include "Config.h"
#include "CrossSection.h"
#include "Absorption.h"
#include "Settings.h"
#include "InputReader.h"
#include "Plotter.h"
#include "Writer.h"

using std::string;
using std::vector;

class Target{
private:
	Settings settings;

	vector<double> z_bins;
	vector<double> energy_boosted;

	// Cross section at rest
	vector< vector<double> > crosssection_at_rest_histogram;

	// Velocity distribution
	vector< vector<double> > velocity_distribution_file;
	vector< vector<double> > velocity_distribution_bins;
	vector< vector<double> > velocity_distribution_histogram;
	vector<double> vdist_norm;
	vector<unsigned int> vdist_centroid;

	// Cross section transformed by velocity distribution
	vector<double> crosssection_histogram;

	// Incident beam
	vector< vector<double> > incident_beam_file;
	vector<double> incident_beam_histogram;

	vector<vector <double> > photon_flux_density_histogram;
	vector<vector <double> > resonance_absorption_density_histogram;

	vector<double> massattenuation_histogram;
	vector<double> transmitted_beam_histogram;

	CrossSection *crossSection;
	Absorption *absorption;
	InputReader *inputReader;
	Plotter *plotter;
	Writer *writer;

	unsigned int target_number;

public:	
	Target(unsigned int ne, unsigned int nz, string name, unsigned int number);
	Target(unsigned int n, Settings &s){ 
		settings = s; 
		target_number = n;
	};
	
	~Target();

	void initialize(vector<double> &energy_bins);
	void calculateCrossSection(vector<double> &energy_bins);
	void calculateIncidentBeam(const vector<double> &energy_bins);
	void calculateIncidentBeam(vector<vector<double> > &photon_flux_density_histogram);

	// Function to modify resonance energies due to Doppler shift
	void boostEnergies();

	// Function to print information to the command line
	void print();

	// Function to plot histograms
	void plot(vector<double> &energy_bins);

	// Functions to calculate histograms
	void calculateCrossSectionAtRest(vector<double> &energy_bins);
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

	// Functions to write histograms to file
	void write(vector<double> &energy_bins);

	// Function to set the incident beam and return the transmitted beam
	double& getTransmittedBeam(){ return transmitted_beam_histogram[0]; };
	void setIncidentBeam(double &trans_beam_bins);

	// Functions to integrate over 2D histograms
	double integrateEZHistogram(vector<double> &energy_bins, vector<double> &z_bins, vector<vector<double> > &ezhist);
	double integrateEEHistogram(vector<double> &energy_bins, vector<vector<double> > &eehist);

	void vDistInfo();
};

#endif
