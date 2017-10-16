#ifndef CROSSSECTION_H 
#define CROSSSECTION_H 1

#include <vector>
#include <string>

#include "math.h"

#include "TCanvas.h"
#include "TLegend.h"

#include "Config.h"

using std::vector;
using std::string;

class CrossSection{
	vector< vector<double> > pconv_crosssection_bins;
	vector< vector<double> > pconv_vdist_bins;

public:
	CrossSection(){};
	
	~CrossSection(){};

	void breit_wigner(vector<double> &energy_bins, vector<double> (&crosssection_bins), double e0, double gamma0, double gamma, double jj, double j0);

	void calculateVelocityBins(vector<double> &energy_bins, vector<double> &velocity_bins, double e0);

	void maxwell_boltzmann(vector<double> &energy_bins, vector<double> &velocity_bins, vector<double> &vdist_bins, vector<double> &params, double mass, double e0);

	void maxwell_boltzmann_approximation(double (&dopplercs_bins)[NBINS], vector<double> &energy_bins, vector< vector<double> > &velocity_bins, vector< vector<double> > &vdist_bins, vector<double> &e0_list, vector<double> &gamma0_list, vector<double> &gamma_list, vector<double> &jj_list, double j0, vector<double> &params, double mass);

	void maxwell_boltzmann_debye(vector<double> &energy_bins, vector<double> &velocity_bins, vector<double> &vdist_bins, vector<double> &params, double mass, double e0);

	void fft_input(vector<double> &energy_bins, vector< vector<double> > &crosssection_bins, vector< vector<double> > &vdist_bins, vector<double> e0_list);

	void integration_input(vector< vector<double> > &crosssection_bins, vector< vector<double> > &vdist_bins);

	void dopplershift(double (&dopplercs_bins)[NBINS], vector<double> &energy_bins, vector< vector <double> > &crosssection_bins, vector< vector<double> > &velocity_bins, vector< vector<double> > &vdist_bins, vector<double> &vdist_norm, vector<double> &e0_list);

	void dopplershiftFFT(double (&dopplercs_bins)[NBINS], vector<double> &energy_bins, vector< vector <double> > &crosssection_bins, vector< vector<double> > &velocity_bins, vector< vector<double> > &vdist_bins, vector<double> &vdist_norm, vector<unsigned int> &vdist_centroid);

	void plot_crosssection(vector<double> &energy_bins, vector< vector<double> > &crosssection, string title, TCanvas* canvas, TLegend* legend, string legend_entry);

	void plot_vdist(vector<double> &velocity_bins, vector<double> &vdist_bins, string title, TCanvas* canvas, TLegend* legend, string legend_entry);

	void plot_dopplershift(vector<double> &energy_bins, vector< vector<double> > &crosssection_bins, double (&dopplercs_bins)[NBINS], string title, TCanvas* canvas, TLegend* legend, string legend_entry);

private:
	double eGamma(double energy, double velocity){
		return (1. + velocity)/(sqrt(1. - velocity*velocity))*energy;
	}

	double delta(double t, double mass){
		return sqrt(2.*kB*t/(mass*AtomicMassUnit));
	}
};

#endif
