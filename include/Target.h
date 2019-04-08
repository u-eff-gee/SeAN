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


#ifndef TARGET_H 
#define TARGET_H 1

#include <string>
#include <vector>

#include "Config.h"
#include "CrossSection.h"
#include "Absorption.h"
#include "Settings.h"
#include "InputReader.h"
#include "Plotter.h"
#include "Writer.h"
#include "Integrator.h"
#include "PhononDensity.h"

using std::string;
using std::vector;

class Target{
private:
	Settings settings;

	vector<double> z_bins;
	vector<double> energy_boosted;

	// Cross section at rest
	vector< vector<double> > crosssection_at_rest_histogram;

	// Velocity distribution
	vector< vector<double> > velocity_distribution_bins;
	vector< vector<double> > velocity_distribution_histogram;
	vector<double> vdist_norm;
	vector<double> velocity_bins_file;
	vector<double> velocity_distribution_file;
	vector<unsigned int> vdist_centroid;

	// Cross section transformed by velocity distribution
	vector<double> energy_bins_file;
	vector<double> cross_section_file;
	vector<double> crosssection_histogram;

	// Phonon density files
	//vector< vector<double> > e_s_file;
	//vector< vector<double> > p_file;
	vector<double> omega_s_file;

	// Incident beam
	vector< vector<double> > incident_beam_file;
	vector<double> incident_beam_histogram;

	vector<vector <double> > photon_flux_density_histogram;
	vector<vector <double> > resonance_absorption_density_histogram;

	vector< vector<double> > mass_attenuation_file;
	vector<double> mass_attenuation_histogram;
	vector<double> transmitted_beam_histogram;

	CrossSection crossSection;
	Absorption absorption;
	InputReader inputReader;
	Plotter plotter;
	Writer writer;
	Integrator integrator;
	PhononDensity phononDensity;

	double n_resonantly_scattered;
	unsigned int target_number;

public:	
	Target(unsigned int n, Settings &s): 
		crossSection(s), 
		absorption(s),
       		inputReader(),
		plotter(s),
		writer(s),
		integrator(),
		phononDensity(s),
		n_resonantly_scattered(0.)
	{ 
		settings = s; 
		target_number = n;
	};
	
	~Target(){};

	// The order of the following functions represents the order in which they are called in a normal calculation. Indented functions are called by the function above.
	void initialize(const vector<double> &energy_bins);
		// Function to modify resonance energies due to Doppler shift
		void boost_and_recoil();
		void calculateCrossSectionAtRest(const vector<double> &energy_bins);
		void calculateVelocityDistribution(const vector<double> &energy_bins);
		void calculateMassAttenuation(const vector<double> &energy_bins);
	void calculateCrossSection(const vector<double> &energy_bins);
	void calculateIncidentBeam(const vector<double> &energy_bins);
	void calculateIncidentBeam(const vector<vector<double> > &photon_flux_density_histogram);
	void calculateTransmission(const vector<double> energy_bins);
	void calculateResonantScattering(const vector<double> energy_bins);

	// Function to print information to the command line
	void print();
	void print_results();

	// Function to plot histograms
	void plot(vector<double> &energy_bins, unsigned int n_setting) ;

	// Functions to calculate histograms

	void calculateTransmittedBeam();
	void calculateZBins();
	void calculateZBins(double z0, double z1);
	void calculatePhotonFluxDensity();
	void calculateResonanceAbsorptionDensity();
	void calculateAbsorption(vector<double> &energy_bins);

	// Functions to write histograms to file
	void write(const vector<double> &energy_bins, const unsigned int n_setting);

	// Functions to write output
	void write_results(string outputfile) const;

	// Function to return the photon flux density. The incident beam on the subsequent target will be read off from this.
	
	vector< vector<double> >& getPhotonFluxDensity(){ return photon_flux_density_histogram; };

	// Functions to integrate over 2D histograms
	double integrateEZHistogram(vector<double> &energy_bins, vector<double> &z_bins, vector<vector<double> > &ezhist);
	double integrateEEHistogram(vector<double> &energy_bins, vector<vector<double> > &eehist);

	void vDistInfo();
};

#endif
