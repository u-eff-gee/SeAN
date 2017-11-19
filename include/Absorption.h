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
#include "Settings.h"

using std::string;
using std::vector;

class Absorption{
private:
	Settings settings;

	vector< vector<double> > matt;

public:
	Absorption(Settings &s){ settings = s; };
	~Absorption(){
		delete &matt;
	};

	// Incident beam calculators
	void const_beam(const vector<double> &energy_bins, vector<double> &incident_beam_histogram);

	void gauss_beam(const vector<double> &energy_bins, vector<double> &incident_beam_histogram);

	void arbitrary_beam(const vector<double> &energy_bins, vector<double> &incident_beam_histogram, const vector< vector<double> > &incident_beam_file);

	// Mass attenuation calculators
	void const_mass_attenuation(vector<double> &mass_attenuation_histogram, const unsigned int target_number);

	void nist_mass_attenuation(const vector<double> &energy_bins, vector<double> &mass_attenuation_histogram, const unsigned int target_number);

	void arbitrary_mass_attenuation(const vector<double> &energy_bins, const vector< vector<double> > mass_attenuation_file, vector<double> &mass_attenuation_histogram);

	void nist_mass_attenuation(const vector<double> &energy_bins, const vector< vector<double> > mass_attenuation_file, vector<double> &mass_attenuation_histogram, const unsigned int target_number);

	// Transmission calculators
	void photon_flux_density(const vector<double> &crosssection_histogram, const vector<double> &mass_attenuation_histogram, const vector<double> &z_bins, const vector<double> &incident_beam_histogram, vector<vector<double> > &photon_flux_density_histogram);

	void resonance_absorption_density(const vector<double> &crosssection_histogram, const vector<vector<double> > &photon_flux_density_histogram, vector< vector<double> > &resonance_absorption_density_histogram);

	void plot_photon_flux_density(vector<double> &energy_bins, vector<double> &z_bins, vector<vector<double> > &photon_flux_density_bins, string title, TCanvas *canvas, TLegend* legend, string legend_entry);

	void plot_resonance_absorption_density(vector<double> &energy_bins, vector<double> &z_bins, vector<vector<double> > &resonance_absorption_density_bins, string title, TCanvas *canvas, TLegend* legend, string legend_entry);

};

#endif
