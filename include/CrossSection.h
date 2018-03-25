#ifndef CROSSSECTION_H 
#define CROSSSECTION_H 1

#include <vector>
#include <string>

#include "math.h"

#include "TCanvas.h"
#include "TLegend.h"

#include "Config.h"
#include "Settings.h"
#include "Integrator.h"

using std::vector;
using std::string;

class CrossSection{

private:
	vector< vector<double> > pconv_crosssection_histogram;
	vector< vector<double> > pconv_velocity_distribution_histogram;

	Settings settings;
	Integrator* integrator;

public:
	CrossSection(Settings &s){ 
		settings = s; 
		integrator = new Integrator();
	};
	
	~CrossSection(){
		delete &pconv_crosssection_histogram;
		delete &pconv_velocity_distribution_histogram;
	};

	// Cross section at rest calculators
	void breit_wigner(const vector<double> &energy_bins, vector< vector<double> > (&crosssection_at_rest_histogram), const vector<double> &energy_boosted, const unsigned int target_number);

	// Velocity distribution calculators
	// Calculator for the velocity bins that correspond to the energy bins
	void calculateVelocityBins(const vector<double> &energy_bins, vector< vector<double> > &velocity_distribution_bins, vector<double> &energy_boosted, unsigned int target_number);

	void absolute_zero(const vector< vector<double> > velocity_distribution_bins, vector< vector<double> > &velocity_distribution_histogram, const unsigned int target_number);

	void maxwell_boltzmann(const vector< vector<double> > &velocity_distribution_bins, vector< vector<double> > &velocity_distribution_histogram, const unsigned int target_number);

	void maxwell_boltzmann_approximation(const vector<double> &energy_bins, vector<double> &crosssection_histogram, const vector<double> &energy_boosted, const unsigned int target_number);
	void maxwell_boltzmann_approximation_debye(const vector<double> &energy_bins, vector<double> &crosssection_histogram, const vector<double> &energy_boosted, const unsigned int target_number);
	const double APPROXIMATION_LIMIT = 0.1;

	double tEff(const double t, const double tD);
	void maxwell_boltzmann_debye(const vector< vector<double> > &velocity_distribution_bins, vector< vector<double> > &velocity_distribution_histogram, const unsigned int target_number);

	void arbitrary_velocity_distribution(const vector< vector<double> > &velocity_distribution_bins, vector< vector<double> > &velocity_distribution_histogram, const vector<double> &velocity_bins_file, const vector<double> &velocity_distribution_file, const vector<double> &energy_boosted, const unsigned int target_number);

	// Cross section calculators
	// Using FFT
	void fft_input(const vector<double> &energy_bins, vector< vector<double> > &crosssection_histogram, vector< vector<double> > &velocity_distribution_histogram, vector<double> energy_boosted);

	void dopplershiftFFT(const vector<double> &energy_bins, vector<double> &crosssection_histogram, vector< vector <double> > &crosssection_at_rest_histogram, vector< vector<double> > &velocity_distribution_bins, vector< vector<double> > &velocity_distribution_histogram, vector<double> &vdist_norm, vector<unsigned int> &vdist_centroid);

	// Using trapezoidal rule
	void integration_input(const vector< vector<double> > &crosssection_bins, const vector< vector<double> > &vdist_bins);

	void dopplershift(const vector<double> &energy_bins, vector<double> &crosssection_histogram, vector< vector <double> > &crosssection_bins, vector< vector<double> > &velocity_distribution_bins, vector< vector<double> > &velocity_distribution_histogram, vector<double> &vdist_norm, vector<double> &energy_boosted);

	// Special cases
	void no_dopplershift(const vector< vector<double> > &crosssection_at_rest_histogram, vector<double> &crosssection_histogram);

	void arbitrary_cross_section(const vector<double> &energy_bins, vector<double> &crosssection_histogram, const vector<double> &energy_bins_file, const vector<double> &cross_section_file);

	// Unit tests
	// Calculate the analytical value integrated cross section, i.e. the integral of d sigma / d E with from 0 to infinity. This is a measure of the correctness of the cross section that is calculated by SeAN
	double integrated_crosssection_analytical(vector<double> energy_boosted, unsigned int target_number) const;

	// Check whether the total integral over the cross section equals the analytical expression for a Breit-Wigner cross section folded with an arbitrary velocity distribution
	void check_crosssection_normalization(const vector<double> &energy_bins, vector<double> &crosssection_histogram, vector<double> energy_boosted, unsigned int target_number) const;

private:
	double eGamma(double energy, double velocity){
		return (1. + velocity)/(sqrt(1. - velocity*velocity))*energy;
	}

	double delta(double t, double mass){
		return sqrt(2.*kB*t/(mass*AtomicMassUnit));
	}
};

class tEff_integrated_function{
	
public:
	double operator()(double t) const {
		return pow(t,3)*(1./(exp(t) - 1.) + 0.5);
	}

	double Eval(double t) const {
		return pow(t,3)*(1./(exp(t) - 1.) + 0.5);
	};
};

#endif
