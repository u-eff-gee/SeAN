#ifndef ABSORPTION_H
#define ABSORPTION_H 1

#include <string>
#include <vector>

#include "Math/Interpolator.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TAxis.h"

#include "Config.h"

using std::string;
using std::vector;

class Absorption{
private:
	vector< vector<double> > matt;
	ROOT::Math::Interpolator *inter;

	unsigned int nbins_matt;
public:
	Absorption(){
		nbins_matt = 0;

		matt.push_back(vector<double>());
		matt.push_back(vector<double>());
	};
	~Absorption(){};

	void read_massattenuation_NIST(double (&energy_bins)[NBINS], double (&massattenuation_bins)[NBINS], string massAttenuation_ID, double mass);

	void absorption(double (&energy_bins)[NBINS], double (&dopplercs_bins)[NBINS]);

	void plot_massattenuation(double (&energy_bins)[NBINS], double (&massattenuation_bins)[NBINS], string title, TCanvas *canvas, TLegend* legend, string legend_entry);

	void plot_total_massattenuation(string title, TCanvas *canvas, TLegend* legend, string legend_entry);
};

#endif
