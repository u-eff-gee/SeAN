/*    
    This file is part of SeAN.

    SeAN is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SeAN is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SeAN.  If not, see <http://www.gnu.org/licenses/>.
*/


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
	~Absorption(){};

	// Incident beam calculators
	void const_beam(vector<double> &incident_beam_histogram);

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

	// Transmission calculators with an additional parameter to be multiplied to the
	// cross section. This is to estimate the uncertainty when the cross section is
	// scaled.
	void photon_flux_density(const vector<double> &crosssection_histogram, const vector<double> &mass_attenuation_histogram, const vector<double> &z_bins, const vector<double> &incident_beam_histogram, vector<vector<double> > &photon_flux_density_histogram, const double cs_enhancement_factor);

	void resonance_absorption_density(const vector<double> &crosssection_histogram, const vector<vector<double> > &photon_flux_density_histogram, vector< vector<double> > &resonance_absorption_density_histogram, const double cs_enhancement_factor);

	// Plot functions
	void plot_photon_flux_density(vector<double> &energy_bins, vector<double> &z_bins, vector<vector<double> > &photon_flux_density_bins, string title, TCanvas *canvas, TLegend* legend, string legend_entry);

	void plot_resonance_absorption_density(vector<double> &energy_bins, vector<double> &z_bins, vector<vector<double> > &resonance_absorption_density_bins, string title, TCanvas *canvas, TLegend* legend, string legend_entry);

};

#endif
