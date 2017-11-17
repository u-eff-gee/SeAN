#include "CrossSection.h"

#include <iostream>
#include <sstream>
#include <algorithm>
#include <complex.h>

#include "omp.h"
#include "fftw3.h"

#include "TGraph.h"
#include "TAxis.h"

using std::cout;
using std::endl;
using std::upper_bound;
using std::max_element;

void CrossSection::breit_wigner(vector<double> &energy_bins, vector<vector<double> > &crosssection_at_rest_bins, vector<double> &energy_boosted, unsigned int target_number){

	for(unsigned int i = 0; i < energy_boosted.size(); ++i){
		double cs_max = PI*0.5*HBARC2/(energy_boosted[i]*energy_boosted[i])*(2.*settings.jj[target_number][i] + 1.)/(2.*settings.ji[target_number] + 1.)*settings.gamma0[target_number][i]*settings.gamma[target_number][i];
		
		for(unsigned int j = 0; j < settings.nbins_e; ++j){
			crosssection_at_rest_bins[i][j] += cs_max / ((energy_bins[j] - energy_boosted[i])*(energy_bins[j] - energy_boosted[i]) + 0.25*settings.gamma[target_number][i]*settings.gamma[target_number][i]);
		}
	}
}

void CrossSection::calculateVelocityBins(vector<double> &energy_bins, vector< vector<double> > &velocity_distribution_bins, vector<double> &energy_boosted, unsigned int target_number){

	// SeAN uses non-equidistant velocity bins. Here, velocity_bins[i] is the velocity that is needed to shift energy_bins[j] to energy_bins[i].
	double energy_ratio = 1.;

	// Find bin of resonance energy
//	long int resonance_position = upper_bound(energy_bins, energy_bins + NBINS, e0) - energy_bins;
//	double resonance_energy = energy_bins[resonance_position];

	for(unsigned int i = 0; i < energy_boosted.size(); ++i){
		for(unsigned int j = 0; j < settings.nbins_e; ++j){
			energy_ratio = energy_bins[j]/energy_boosted[i];
			velocity_distribution_bins[i][j] = (1. - energy_ratio*energy_ratio)/(1. + energy_ratio*energy_ratio);
		}
	}
}

void CrossSection::maxwell_boltzmann(vector<double> &energy_bins, vector< vector<double> > &velocity_distribution_bins, vector< vector<double> > &velocity_distribution_histogram, vector<double> &energy_boosted, unsigned int target_number){

	double c1 = sqrt(settings.mass[target_number]*AtomicMassUnit/(2*PI*kB*settings.vDistParams[target_number][0]));
      	double c2 = -1./pow(delta(settings.vDistParams[target_number][0], settings.mass[target_number]), 2);

      	for(unsigned int i = 0; i < energy_boosted.size(); ++i){
		for(unsigned int j = 0; j < settings.nbins_e; ++j){
			velocity_distribution_histogram[i][j] = c1*exp(c2*velocity_distribution_bins[i][j]*velocity_distribution_bins[i][j]);
		}
	}
}

void CrossSection::maxwell_boltzmann_debye(vector<double> &energy_bins, vector<double> &velocity_bins, vector<double> &vdist_bins, vector<double> &params, double mass, double e0){
	;
}

void CrossSection::integration_input(vector< vector<double> > &crosssection_histogram, vector< vector<double> > &vdist_bins){
	// Zero-pad the cross section array with nbins_e bins on both sides, to get an array of size 3*nbins_e. This is necessary, because the integration works by shifting the two arrays against each other.
	for(unsigned int i = 0; i < crosssection_histogram.size(); ++i){
		pconv_crosssection_histogram.push_back(vector<double>(3*settings.nbins_e, 0.));

		for(unsigned int j = 0; j < settings.nbins_e; ++j){
			pconv_crosssection_histogram[i][j + settings.nbins_e] = crosssection_histogram[i][j];
		}
	}
}

void CrossSection::dopplershift(vector<double> &energy_bins, vector<double> &crosssection_histogram, vector< vector <double> > &crosssection_bins, vector< vector<double> > &velocity_distribution_bins, vector< vector<double> > &velocity_distribution_histogram, vector<double> &vdist_norm, vector<double> &energy_boosted){

	// Integrate using trapezoidal rule
	
	double resonance_energy_squared = 1.;
	
	for(unsigned int n = 0; n < pconv_crosssection_histogram.size(); ++n){
		long unsigned int resonance_position = (unsigned long int) (upper_bound(energy_bins.begin(), energy_bins.end(), energy_boosted[n]) - energy_bins.begin());
		resonance_energy_squared = energy_boosted[n]*energy_boosted[n];

		for(unsigned int i = 0; i < settings.nbins_e; ++i){
			for(unsigned int j = 0; j < settings.nbins_e - 1; ++j){
				crosssection_histogram[i] += 0.5*(
						velocity_distribution_histogram[n][j]*pconv_crosssection_histogram[n][settings.nbins_e + i + resonance_position - j]*resonance_energy_squared/(energy_bins[j]*energy_bins[j]) 
						+ velocity_distribution_histogram[n][j + 1]*pconv_crosssection_histogram[n][settings.nbins_e + i + resonance_position - (j + 1)]*resonance_energy_squared/(energy_bins[j + 1]*energy_bins[j + 1]))
					*(velocity_distribution_bins[n][j] - velocity_distribution_bins[n][j + 1]); 
			}
		}
	}
}

void CrossSection::fft_input(vector<double> &energy_bins, vector< vector<double> > &crosssection_at_rest_histogram, vector< vector<double> > &velocity_distribution_histogram, vector<double> energy_boosted){
	// Preparatory work for the FFT convolution:
	// Convert the velocity distribution w(v) to w(E)*dv/dE so that both the cross section and w depend on the energy and the integral can be calculated
	// Multiply the infinitesimal integration interval dE to the cross section bin
	// Store the converted cross section and velocity distribution histograms in new vectors, since fftw destroys its input arrays on execution
	double de = energy_bins[1] - energy_bins[0];
	double energy_ratio_squared = 1.;

	for(unsigned int i = 0; i < energy_boosted.size(); ++i){
		pconv_crosssection_histogram.push_back(vector<double>(settings.nbins_e, 0.));
		pconv_velocity_distribution_histogram.push_back(vector<double>(settings.nbins_e, 0.));

		for(unsigned int j = 0; j < settings.nbins_e; ++j){
			pconv_crosssection_histogram[i][j] = crosssection_at_rest_histogram[i][j]*de;

			energy_ratio_squared = energy_bins[j]*energy_bins[j]/(energy_boosted[i]*energy_boosted[i]);	
			pconv_velocity_distribution_histogram[i][j] = -velocity_distribution_histogram[i][j]*(-4.*energy_ratio_squared/energy_boosted[i])/((1. + energy_ratio_squared)*(1. + energy_ratio_squared));
		}
	}
}

void CrossSection::dopplershiftFFT(vector<double> &energy_bins, vector<double> &crosssection_histogram, vector< vector <double> > &crosssection_at_rest_histogram, vector< vector<double> > &velocity_distribution_bins, vector< vector<double> > &velocity_distribution_histogram, vector<double> &vdist_norm, vector<unsigned int> &vdist_centroid){

	fftw_plan vdist_plan, crosssection_plan, product_fft_plan;
	fftw_complex vdist_fft[settings.nbins_e/2 + 1] = {{0.}};
	fftw_complex crosssection_fft[settings.nbins_e/2 + 1] = {{0.}};
	fftw_complex product_fft[settings.nbins_e/2 + 1] = {{0.}};
	double convolution[settings.nbins_e] = {0.};

	for(unsigned int i = 0; i < crosssection_at_rest_histogram.size(); ++i){

		vdist_plan = fftw_plan_dft_r2c_1d((int) settings.nbins_e, &pconv_velocity_distribution_histogram[i][0], vdist_fft, FFTW_ESTIMATE);
		crosssection_plan = fftw_plan_dft_r2c_1d((int) settings.nbins_e, &pconv_crosssection_histogram[i][0], crosssection_fft, FFTW_ESTIMATE);

		// Calculate Fourier transform of both arrays
		fftw_execute(vdist_plan);
		fftw_execute(crosssection_plan);

		// Multiply the (complex) Fourier - transformed arrays
		for(unsigned int j = 0; j < settings.nbins_e / 2 + 1; ++j){
			product_fft[j][0] = crosssection_fft[j][0]*vdist_fft[j][0] - crosssection_fft[j][1]*vdist_fft[j][1];
			product_fft[j][1] = crosssection_fft[j][0]*vdist_fft[j][1] + crosssection_fft[j][1]*vdist_fft[j][0];
		}

		// Transform the product back
		product_fft_plan = fftw_plan_dft_c2r_1d((int) settings.nbins_e, product_fft, convolution, FFTW_ESTIMATE);
		fftw_execute(product_fft_plan);

		 // Add the i-th convoluted cross section to the total cross section, observing that the transformation back changed the order of the elements in the array and leaves them scaled by settings.nbins_e.
		unsigned int shift = vdist_centroid[i];
		double inverse_norm = (double) 1./(vdist_norm[i]*settings.nbins_e);
		shift = settings.nbins_e - shift;
		for(unsigned int j = 0; j < settings.nbins_e; ++j){
			if(j + shift < settings.nbins_e){
				crosssection_histogram[j + shift] += convolution[j]*inverse_norm;
			} else{
				crosssection_histogram[j - (settings.nbins_e - shift)] += convolution[j]*inverse_norm;
			}
		}
	}
}

void CrossSection::maxwell_boltzmann_approximation(vector<double> &dopplercs_bins, vector<double> &energy_bins, vector< vector<double> > &velocity_bins, vector< vector<double> > &vdist_bins, vector<double> &e0_list, vector<double> &gamma0_list, vector<double> &gamma_list, vector<double> &jj_list, double j0, vector<double> &params, double mass){

// Calculate velocity distribution as in maxwell_boltzmann

	for(unsigned int i = 0; i < e0_list.size(); ++i){
		double c1 = sqrt(mass*AtomicMassUnit/(2*PI*kB*params[0]));
		double c2 = -1./pow(delta(params[0], mass), 2);

		for(unsigned int j = 0; j < NBINS; ++j){
			vdist_bins[i][j] = c1*exp(c2*velocity_bins[i][j]*velocity_bins[i][j]);
		}
	}

// Calculate doppler-shifted cross section directly

	for(unsigned int i = 0; i < e0_list.size(); ++i){
		double doppler_width = sqrt(2.*kB*params[0]/(mass*AtomicMassUnit))*e0_list[i];

		if(gamma_list[i]/doppler_width > APPROXIMATION_LIMIT){
			cout << "> Warning: CrossSection.cpp: maxwell_approximation(): Gamma/Delta = " << gamma_list[i]/doppler_width << " > " << APPROXIMATION_LIMIT << ", the applicability of the approximation for the doppler-shifted cross section may not be good." << endl;
			cout << "\tE0 = " << e0_list[i] << " eV" << endl;
			cout << "\tGAMMA = " << gamma_list[i] << " eV" << endl;
			cout << "\tMASS = " << mass << " u" << endl;
			cout << "\tTEFF= " << params[0] << " K" << endl;
		}

		double cs_max = 2.*PI*HBARC2/(e0_list[i]*e0_list[i])*(2.*jj_list[i] + 1.)/(2. * j0 + 1.)*gamma0_list[i]/gamma_list[i]*sqrt(PI)/(2.*doppler_width/gamma_list[i]);
		
		for(unsigned int j = 0; j < NBINS; ++j){
			dopplercs_bins[j] += cs_max*exp(-(energy_bins[j] - e0_list[i])*(energy_bins[j] - e0_list[i])/(doppler_width*doppler_width));
		}
	}
}
