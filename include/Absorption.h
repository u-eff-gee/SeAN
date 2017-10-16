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

	void read_massattenuation_NIST(vector<double> &energy_bins, vector<double> &massattenuation_bins, string massAttenuation_ID, double mass);

	void const_beam(vector<double> &energy_bins, vector<double> &incident_beam_bins, vector<double> beamParams);

	void gauss_beam(vector<double> &energy_bins, vector<double> &incident_beam_bins, vector<double> beamParams);

	void photon_flux_density(vector<double> &dopplercs_bins, vector<double> &massattenuation_bins, vector<double> &z_bins, vector<double> &incident_beam_bins, vector<vector<double> > &photon_flux_density_bins);

	void resonance_absorption_density(vector<double> &dopplercs_bins, vector<vector<double> > &photon_flux_density_bins, vector<vector<double> > &resonance_absorption_density_bins);

	void plot_massattenuation(vector<double> &energy_bins, vector<double> &massattenuation_bins, string title, TCanvas *canvas, TLegend* legend, string legend_entry);

	void plot_total_massattenuation(string title, TCanvas *canvas, TLegend* legend, string legend_entry);

	void plot_photon_flux_density(vector<double> &energy_bins, vector<double> &z_bins, vector<vector<double> > &photon_flux_density_bins, string title, TCanvas *canvas, TLegend* legend, string legend_entry);

	void plot_resonance_absorption_density(vector<double> &energy_bins, vector<double> &z_bins, vector<vector<double> > &resonance_absorption_density_bins, string title, TCanvas *canvas, TLegend* legend, string legend_entry);

};

#endif
