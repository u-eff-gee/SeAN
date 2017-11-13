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
	vector< vector<double> > pconv_crosssection_bins;
	vector< vector<double> > pconv_vdist_bins;

	Settings settings;

public:
	CrossSection(Settings &s){ settings = s; };
	
	~CrossSection(){
		delete &pconv_crosssection_bins;
		delete &pconv_vdist_bins;
	};

	// Cross section at rest calculators
	void breit_wigner(vector<double> &energy_bins, vector< vector<double> > (&crosssection_at_rest_histogram), vector<double> &energy_boosted, unsigned int target_number);

	// Velocity distribution calculators
	// Calculator for the velocity bins that correspond to the energy bins
	void calculateVelocityBins(vector<double> &energy_bins, vector< vector<double> > &velocity_distribution_bins, vector<double> &energy_boosted, unsigned int target_number);

	void maxwell_boltzmann(vector<double> &energy_bins, vector< vector<double> > &velocity_distribution_bins, vector< vector<double> > &velocity_distribution_histogram, vector<double> &energy_boosted, unsigned int target_number);

	void maxwell_boltzmann_approximation(vector<double> &dopplercs_bins, vector<double> &energy_bins, vector< vector<double> > &velocity_bins, vector< vector<double> > &vdist_bins, vector<double> &e0_list, vector<double> &gamma0_list, vector<double> &gamma_list, vector<double> &jj_list, double j0, vector<double> &params, double mass);

	void maxwell_boltzmann_debye(vector<double> &energy_bins, vector<double> &velocity_bins, vector<double> &vdist_bins, vector<double> &params, double mass, double e0);

	void fft_input(vector<double> &energy_bins, vector< vector<double> > &crosssection_bins, vector< vector<double> > &vdist_bins, vector<double> e0_list);

	void integration_input(vector< vector<double> > &crosssection_bins, vector< vector<double> > &vdist_bins);

	void dopplershift(vector<double> &dopplercs_bins, vector<double> &energy_bins, vector< vector <double> > &crosssection_bins, vector< vector<double> > &velocity_bins, vector< vector<double> > &vdist_bins, vector<double> &vdist_norm, vector<double> &e0_list);

	void dopplershiftFFT(vector<double> &dopplercs_bins, vector<double> &energy_bins, vector< vector <double> > &crosssection_bins, vector< vector<double> > &velocity_bins, vector< vector<double> > &vdist_bins, vector<double> &vdist_norm, vector<unsigned int> &vdist_centroid);

	void plot_crosssection(vector<double> &energy_bins, vector< vector<double> > &crosssection, string title, TCanvas* canvas, TLegend* legend, string legend_entry);

	void plot_vdist(vector<double> &velocity_bins, vector<double> &vdist_bins, string title, TCanvas* canvas, TLegend* legend, string legend_entry);

	void plot_dopplershift(vector<double> &energy_bins, vector< vector<double> > &crosssection_bins, vector<double> &dopplercs_bins, string title, TCanvas* canvas, TLegend* legend, string legend_entry);

private:
	double eGamma(double energy, double velocity){
		return (1. + velocity)/(sqrt(1. - velocity*velocity))*energy;
	}

	double delta(double t, double mass){
		return sqrt(2.*kB*t/(mass*AtomicMassUnit));
	}
};

#endif
