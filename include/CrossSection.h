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

	void breit_wigner(double (&energy_bins)[NBINS_E], double (&crosssection_bins)[NBINS_V], vector<double> &e0_list, vector<double> &gamma0_list, vector<double> &gamma_list, vector<double> &jj_list, double j0);

	void maxwell_boltzmann(double (&velocity_bins)[NBINS_V], double (&vdist_bins)[NBINS_V], vector<double> &params, double mass);

	void maxwell_boltzmann_debye(double (&velocity_bins)[NBINS_V], double (&vdist_bins)[NBINS_E], vector<double> &params, double mass);

	void dopplershift(double (&dopplercs_bins)[NBINS_E], double (&energy_bins)[NBINS_E], double (&crosssection_bins)[NBINS_V], double (&velocity_bins)[NBINS_V], double (&vdist_bins)[NBINS_E]);

	void plot_crosssection(double (&energies)[NBINS_E], double (&crosssection)[NBINS_E], string title, TCanvas* canvas, TLegend* legend, string legend_entry, bool add_to_existing_canvas);

	void plot_vdist(double (&energies)[NBINS_E], double (&crosssection)[NBINS_E], string title, TCanvas* canvas, TLegend* legend, string legend_entry, bool add_to_existing_canvas);

private:
	double eGamma(double energy, double velocity){
		return (1. + velocity)/(sqrt(1. - velocity*velocity))*energy;
	}

	double delta(double t, double mass){
		return sqrt(kB*t/(mass*AtomicMassUnit));
	}
};

#endif
