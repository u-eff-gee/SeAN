#ifndef CROSSSECTION_H 
#define CROSSSECTION_H 1

#include <vector>
#include <string>

#include "math.h"

#include "TCanvas.h"
#include "TLegend.h"

#include "Config.h"
#include "Settings.h"

using std::vector;
using std::string;

class CrossSection{

private:
	vector< vector<double> > pconv_crosssection_histogram;
	vector< vector<double> > pconv_velocity_distribution_histogram;

	Settings settings;

public:
	CrossSection(Settings &s){ settings = s; };
	
	~CrossSection(){
		delete &pconv_crosssection_histogram;
		delete &pconv_velocity_distribution_histogram;
	};

	// Cross section at rest calculators
	void breit_wigner(vector<double> &energy_bins, vector< vector<double> > (&crosssection_at_rest_histogram), vector<double> &energy_boosted, unsigned int target_number);

	// Velocity distribution calculators
	// Calculator for the velocity bins that correspond to the energy bins
	void calculateVelocityBins(vector<double> &energy_bins, vector< vector<double> > &velocity_distribution_bins, vector<double> &energy_boosted, unsigned int target_number);

	void maxwell_boltzmann(vector< vector<double> > &velocity_distribution_bins, vector< vector<double> > &velocity_distribution_histogram, unsigned int target_number);

	void maxwell_boltzmann_approximation(vector<double> &energy_bins, vector<double> &dopplercs_bins, vector< vector<double> > &velocity_bins, vector< vector<double> > &vdist_bins, vector<double> &e0_list, vector<double> &gamma0_list, vector<double> &gamma_list, vector<double> &jj_list, double j0, vector<double> &params, double mass);

	void maxwell_boltzmann_debye(vector<double> &velocity_bins, vector<double> &vdist_bins, vector<double> &params, double mass, double e0);

	void arbitrary_velocity_distribution(vector< vector<double> > &velocity_distribution_bins, vector< vector<double> > &velocity_distribution_histogram, vector< vector<double> > &velocity_distribution_file, vector<double> &energy_boosted, unsigned int target_number);

	// Cross section calculators
	// Using FFT
	void fft_input(vector<double> &energy_bins, vector< vector<double> > &crosssection_histogram, vector< vector<double> > &velocity_distribution_histogram, vector<double> energy_boosted);

	void dopplershiftFFT(vector<double> &energy_bins, vector<double> &crosssection_histogram, vector< vector <double> > &crosssection_at_rest_histogram, vector< vector<double> > &velocity_distribution_bins, vector< vector<double> > &velocity_distribution_histogram, vector<double> &vdist_norm, vector<unsigned int> &vdist_centroid);

	// Using trapezoidal rule
	void integration_input(vector< vector<double> > &crosssection_bins, vector< vector<double> > &vdist_bins);

	void dopplershift(vector<double> &energy_bins, vector<double> &crosssection_histogram, vector< vector <double> > &crosssection_bins, vector< vector<double> > &velocity_distribution_bins, vector< vector<double> > &velocity_distribution_histogram, vector<double> &vdist_norm, vector<double> &energy_boosted);

private:
	double eGamma(double energy, double velocity){
		return (1. + velocity)/(sqrt(1. - velocity*velocity))*energy;
	}

	double delta(double t, double mass){
		return sqrt(2.*kB*t/(mass*AtomicMassUnit));
	}
};

#endif
