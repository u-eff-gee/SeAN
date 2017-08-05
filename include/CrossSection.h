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
public:
	CrossSection(){};
	
	~CrossSection(){};

	void breit_wigner(double (&energy_bins)[NBINS], double (&crosssection_bins)[NBINS], vector<double> &e0_list, vector<double> &gamma0_list, vector<double> &gamma_list, vector<double> &jj_list, double j0);

	void maxwell_boltzmann(double (&energy_bins)[NBINS], vector<double> &velocity_bins, vector<double> &vdist_bins, vector<double> &params, double mass, double e0);

	void maxwell_boltzmann_debye(double (&energy_bins)[NBINS], vector<double> &velocity_bins, vector<double> &vdist_bins, vector<double> &params, double mass, double e0);

	void dopplershift(double (&dopplercs_bins)[NBINS], double (&energy_bins)[NBINS], double (&crosssection_bins)[NBINS], vector<double> &velocity_bins, vector<double> &vdist_bins);

	void plot_crosssection(double (&energies)[NBINS], double (&crosssection)[NBINS], string title, TCanvas* canvas, TLegend* legend, string legend_entry, bool add_to_existing_canvas);

	void plot_vdist(vector<double> &velocity_bins, vector<double> &vdist_bins, string title, TCanvas* canvas, TLegend* legend, string legend_entry, bool add_to_existing_canvas);

private:
	double eGamma(double energy, double velocity){
		return (1. + velocity)/(sqrt(1. - velocity*velocity))*energy;
	}

	double delta(double t, double mass){
		return sqrt(kB*t/(mass*AtomicMassUnit));
	}

	// Binary search to find the bin that corresponds to a given floating-point value
	unsigned int get_bin(double d, double (&array)[NBINS]){
		int start = 0;
		int stop = NBINS - 1;

		while(start != stop){
			if(array[stop/2] < d){
				stop = start;
			} else{
				start = stop/2;
			}
		}

		return start;
	}
};

#endif
