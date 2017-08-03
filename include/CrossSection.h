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

	void breit_wigner(double (&energies)[NBINS_E], double (&crosssection)[NBINS_E], vector<double> &e0, vector<double> &gamma0, vector<double> &gamma, vector<double> &jj, double j0);

	void absolute_zero(){};

	void maxwell_boltzmann(double (&energies)[NBINS_E], double (&crosssection)[NBINS_E], vector<double> params);

	void maxwell_boltzmann_debye(double (&energies)[NBINS_E], double (&crosssection)[NBINS_E], vector<double> params);

	void plot(double (&energies)[NBINS_E], double (&crosssection)[NBINS_E], string title, TCanvas* canvas, TLegend* legend, string legend_entry, bool add_to_existing_canvas);

private:
	double eGamma(double energy, double velocity){
		return (1. + velocity)/(sqrt(1. - velocity*velocity))*energy;
	}
};

#endif
