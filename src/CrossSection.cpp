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


#include "CrossSection.h"

#include <iostream>
#include <sstream>
#include <algorithm>

#include "omp.h"
#include "fftw3.h"

#include "Math/Interpolator.h"
#include "Math/Integrator.h"
#include "Math/Functor.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"

using std::cout;
using std::endl;
using std::upper_bound;
using std::max_element;

void CrossSection::breit_wigner(const vector<double> &energy_bins, vector<vector<double> > &crosssection_at_rest_bins, const vector<double> &energy_boosted, const unsigned int target_number){

	for(unsigned int i = 0; i < energy_boosted.size(); ++i){
		double cs_max = PI*0.5*HBARC2/(energy_boosted[i]*energy_boosted[i])*(2.*settings.jj[target_number][i] + 1.)/(2.*settings.ji[target_number] + 1.)*settings.gamma0[target_number][i]*settings.gamma[target_number][i];
		
		for(unsigned int j = 0; j < settings.nbins_e; ++j){
			crosssection_at_rest_bins[i][j] += cs_max / ((energy_bins[j] - energy_boosted[i])*(energy_bins[j] - energy_boosted[i]) + 0.25*settings.gamma[target_number][i]*settings.gamma[target_number][i]);
		}
	}
}

void CrossSection::calculateVelocityBins(const vector<double> &energy_bins, vector< vector<double> > &velocity_distribution_bins, vector<double> &energy_boosted, unsigned int target_number){

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

void CrossSection::absolute_zero(const vector< vector<double> > velocity_distribution_bins, vector< vector<double> > &velocity_distribution_histogram, const unsigned int target_number){
	// At absolute zero (in the classical sense), the nucleus is completely at rest.
	// I.e., the velocity distribution represents a Dirac delta function.
	// Find the velocity distribution bin that corresponds to zero, and fill only that histogram bin.
	
	// Inefficient search algorithm! Use std::algorithm::upper_bound or similar algorithm!

	unsigned int zero_bin = (unsigned int) settings.nbins_e/2;
	double absolute_value = 0.;

	for(unsigned int i = 0; i < velocity_distribution_bins.size(); ++i){
		absolute_value = fabs(velocity_distribution_bins[i][0]);
		zero_bin = 0;
		for(unsigned int j = 0; j < velocity_distribution_bins[i].size(); ++j){
			if(fabs(velocity_distribution_bins[i][j]) < absolute_value){
				absolute_value = fabs(velocity_distribution_bins[i][j]);
				zero_bin = j;
			}
		}

		velocity_distribution_histogram[i][zero_bin] = 1./fabs(velocity_distribution_bins[i][zero_bin + 1] - velocity_distribution_bins[i][zero_bin]);
	}
}

void CrossSection::maxwell_boltzmann(const vector< vector<double> > &velocity_distribution_bins, vector< vector<double> > &velocity_distribution_histogram, const unsigned int target_number){

	double c1 = sqrt(settings.mass[target_number]*AtomicMassUnit/(2*PI*kB*settings.dopplerParams[target_number][0]));
      	double c2 = -1./pow(delta(settings.dopplerParams[target_number][0], settings.mass[target_number]), 2);

      	for(unsigned int i = 0; i < velocity_distribution_bins.size(); ++i){
		for(unsigned int j = 0; j < settings.nbins_e; ++j){
			velocity_distribution_histogram[i][j] = c1*exp(c2*velocity_distribution_bins[i][j]*velocity_distribution_bins[i][j]);
		}
	}
}

double CrossSection::tEff(const double t, const double tD){
	
	tEff_integrated_function f;
	ROOT::Math::Functor1D wf(f);
	ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE);
	ig.SetFunction(wf);
	
	return 3.*pow(t, 4)/pow(tD, 3)*ig.Integral(0, tD/t);
}

void CrossSection::maxwell_boltzmann_debye(const vector< vector<double> > &velocity_distribution_bins, vector< vector<double> > &velocity_distribution_histogram, const unsigned int target_number){

// Calculate effective temperature in debye approximation
	double teff = tEff(settings.dopplerParams[target_number][0], settings.dopplerParams[target_number][1]);

	double c1 = sqrt(settings.mass[target_number]*AtomicMassUnit/(2*PI*kB*teff));
      	double c2 = -1./pow(delta(teff, settings.mass[target_number]), 2);

      	for(unsigned int i = 0; i < velocity_distribution_bins.size(); ++i){
		for(unsigned int j = 0; j < settings.nbins_e; ++j){
			velocity_distribution_histogram[i][j] = c1*exp(c2*velocity_distribution_bins[i][j]*velocity_distribution_bins[i][j]);
		}
	}
}

void CrossSection::arbitrary_velocity_distribution(const vector< vector<double> > &velocity_distribution_bins, vector< vector<double> > &velocity_distribution_histogram, const vector<double> &velocity_bins_file, const vector<double> &velocity_distribution_file, const vector<double> &energy_boosted, const unsigned int target_number){
	// Interpolate data from file	
	ROOT::Math::Interpolator inter((unsigned int) velocity_bins_file.size(), ROOT::Math::Interpolation::kCSPLINE);
	inter.SetData((unsigned int) velocity_bins_file.size(), &velocity_bins_file[0], &velocity_distribution_file[0]);

	// Evaluate at the given velocity bins
	for(unsigned int i = 0; i < velocity_distribution_bins.size(); ++i){
		for(unsigned int j = 0; j < settings.nbins_e; ++j){
			velocity_distribution_histogram[i][j] = inter.Eval(velocity_distribution_bins[i][j]);		
		}
	}
}

void CrossSection::arbitrary_cross_section(const vector<double> &energy_bins, vector<double> &crosssection_histogram, const vector<double> &energy_bins_file, const vector<double> &cross_section_file){

	// Interpolate data from file	
	ROOT::Math::Interpolator inter((unsigned int) energy_bins_file.size(), ROOT::Math::Interpolation::kCSPLINE);
	inter.SetData((unsigned int) energy_bins_file.size(), &energy_bins_file[0], &cross_section_file[0]);

	// Evaluate at the given cross section bins
	for(unsigned int i = 0; i < settings.nbins_e; ++i){
			crosssection_histogram[i] = inter.Eval(energy_bins[i]);		
	}
}

void CrossSection::integration_input(const vector< vector<double> > &crosssection_histogram, const vector< vector<double> > &vdist_bins){
	// Zero-pad the cross section array with nbins_e bins on both sides, to get an array of size 3*nbins_e. This is necessary, because the integration works by shifting the two arrays against each other.
	for(unsigned int i = 0; i < crosssection_histogram.size(); ++i){
		pconv_crosssection_histogram.push_back(vector<double>(3*settings.nbins_e, 0.));

		for(unsigned int j = 0; j < settings.nbins_e; ++j){
			pconv_crosssection_histogram[i][j + settings.nbins_e] = crosssection_histogram[i][j];
		}
	}
}

void CrossSection::dopplershift(const vector<double> &energy_bins, vector<double> &crosssection_histogram, vector< vector <double> > &crosssection_bins, vector< vector<double> > &velocity_distribution_bins, vector< vector<double> > &velocity_distribution_histogram, vector<double> &vdist_norm, vector<double> &energy_boosted){

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

void CrossSection::no_dopplershift(const vector< vector<double> > &crosssection_at_rest_histogram, vector<double> &crosssection_histogram){

	// Simply sum up the cross sections of the single resonances without changing their shape
	for(unsigned int i = 0; i < crosssection_at_rest_histogram.size(); ++i){
		for(unsigned int j = 0; j < crosssection_at_rest_histogram[i].size(); ++j){
			crosssection_histogram[j] += crosssection_at_rest_histogram[i][j];
		}
	}
}

void CrossSection::fft_input(const vector<double> &energy_bins, vector< vector<double> > &crosssection_at_rest_histogram, vector< vector<double> > &velocity_distribution_histogram, vector<double> energy_boosted){
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

void CrossSection::dopplershiftFFT(const vector<double> &energy_bins, vector<double> &crosssection_histogram, vector< vector <double> > &crosssection_at_rest_histogram, vector< vector<double> > &velocity_distribution_bins, vector< vector<double> > &velocity_distribution_histogram, vector<double> &vdist_norm, vector<unsigned int> &vdist_centroid){

	fftw_plan vdist_plan, crosssection_plan, product_fft_plan;
	vector<fftw_complex> vdist_fft(settings.nbins_e/2 + 1);
	vector<fftw_complex> crosssection_fft(settings.nbins_e/2 + 1);
	vector<fftw_complex> product_fft(settings.nbins_e/2 + 1);
	vector<double> convolution(settings.nbins_e, 0.);

	for(unsigned int i = 0; i < crosssection_at_rest_histogram.size(); ++i){

		vdist_plan = fftw_plan_dft_r2c_1d((int) settings.nbins_e, &pconv_velocity_distribution_histogram[i][0], &vdist_fft[0], FFTW_ESTIMATE);
		crosssection_plan = fftw_plan_dft_r2c_1d((int) settings.nbins_e, &pconv_crosssection_histogram[i][0], &crosssection_fft[0], FFTW_ESTIMATE);

		// Calculate Fourier transform of both arrays
		fftw_execute(vdist_plan);
		fftw_execute(crosssection_plan);

		// Multiply the (complex) Fourier - transformed arrays
		for(unsigned int j = 0; j < settings.nbins_e / 2 + 1; ++j){
			product_fft[j][0] = crosssection_fft[j][0]*vdist_fft[j][0] - crosssection_fft[j][1]*vdist_fft[j][1];
			product_fft[j][1] = crosssection_fft[j][0]*vdist_fft[j][1] + crosssection_fft[j][1]*vdist_fft[j][0];
		}

		// Transform the product back
		product_fft_plan = fftw_plan_dft_c2r_1d((int) settings.nbins_e, &product_fft[0], &convolution[0], FFTW_ESTIMATE);
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

void CrossSection::maxwell_boltzmann_approximation(const vector<double> &energy_bins, vector<double> &crosssection_histogram, const vector<double> &energy_boosted, const unsigned int target_number){

 // Calculate doppler-shifted cross section directly
	for(unsigned int i = 0; i < energy_boosted.size(); ++i){
		double doppler_width = sqrt(2.*kB*settings.dopplerParams[target_number][0]/(settings.mass[target_number]*AtomicMassUnit))*energy_boosted[i];

		if(settings.gamma[target_number][i]/doppler_width > APPROXIMATION_LIMIT){
			cout << "Warning: " << __FILE__ << ":" << __LINE__ << ": "; 
			cout << "maxwell_boltzmann_approximation(): Gamma/Delta = " << settings.gamma[target_number][i]/doppler_width << " > " << APPROXIMATION_LIMIT << ", the approximation of the doppler-shifted cross section may not be good." << endl;
			cout << "\tE0 = " << energy_boosted[i] << " eV" << endl;
			cout << "\tGAMMA = " << settings.gamma[target_number][i] << " eV" << endl;
			cout << "\tMASS = " << settings.mass[target_number] << " u" << endl;
			cout << "\tTEFF= " << settings.dopplerParams[target_number][0] << " K" << endl;
		}

		double cs_max = 2.*PI*HBARC2/(energy_boosted[i]*energy_boosted[i])*(2.*settings.jj[target_number][i] + 1.)/(2. * settings.ji[target_number] + 1.)*settings.gamma0[target_number][i]/settings.gamma[target_number][i]*sqrt(PI)/(2.*doppler_width/settings.gamma[target_number][i]);
		
		for(unsigned int j = 0; j < settings.nbins_e; ++j){
			crosssection_histogram[j] += cs_max*exp(-(energy_bins[j] - energy_boosted[i])*(energy_bins[j] - energy_boosted[i])/(doppler_width*doppler_width));
		}
	}
}

void CrossSection::maxwell_boltzmann_approximation_debye(const vector<double> &energy_bins, vector<double> &crosssection_histogram, const vector<double> &energy_boosted, const unsigned int target_number){

// Calculate effective temperature in debye approximation
	double teff = tEff(settings.dopplerParams[target_number][0], settings.dopplerParams[target_number][1]);

 // Calculate doppler-shifted cross section directly
	for(unsigned int i = 0; i < energy_boosted.size(); ++i){
		double doppler_width = sqrt(2.*kB*teff/(settings.mass[target_number]*AtomicMassUnit))*energy_boosted[i];

		if(settings.gamma[target_number][i]/doppler_width > APPROXIMATION_LIMIT){
			cout << "Warning: " << __FILE__ << ":" << __LINE__ << ": "; 
			cout << "maxwell_boltzmann_approximation(): Gamma/Delta = " << settings.gamma[target_number][i]/doppler_width << " > " << APPROXIMATION_LIMIT << ", the approximation of the doppler-shifted cross section may not be good." << endl;
			cout << "\tE0 = " << energy_boosted[i] << " eV" << endl;
			cout << "\tGAMMA = " << settings.gamma[target_number][i] << " eV" << endl;
			cout << "\tMASS = " << settings.mass[target_number] << " u" << endl;
			cout << "\tTEFF= " << teff << " K" << endl;
		}

		double cs_max = 2.*PI*HBARC2/(energy_boosted[i]*energy_boosted[i])*(2.*settings.jj[target_number][i] + 1.)/(2. * settings.ji[target_number] + 1.)*settings.gamma0[target_number][i]/settings.gamma[target_number][i]*sqrt(PI)/(2.*doppler_width/settings.gamma[target_number][i]);
		
		for(unsigned int j = 0; j < settings.nbins_e; ++j){
			crosssection_histogram[j] += cs_max*exp(-(energy_bins[j] - energy_boosted[i])*(energy_bins[j] - energy_boosted[i])/(doppler_width*doppler_width));
		}
	}
}

double CrossSection::integrated_crosssection_analytical(vector<double> energy_boosted, unsigned int target_number) const {
	double cs_int = 0.;

	for(unsigned i = 0; i < energy_boosted.size(); ++i){
		cs_int += PI2*HBARC2/(energy_boosted[i]*energy_boosted[i])*(2.*settings.jj[target_number][i] + 1.)/(2.*settings.ji[target_number] + 1.)*settings.gamma0[target_number][i];
	}

	return cs_int;
}
	
void CrossSection::check_crosssection_normalization(const vector<double> &energy_bins, vector<double> &crosssection_histogram, vector<double> energy_boosted, unsigned int target_number) const {

	double cs_int_analytical = integrated_crosssection_analytical(energy_boosted, target_number);
	double cs_int_numerical = integrator->trapezoidal_rule(energy_bins, crosssection_histogram);

	cout << HORIZONTAL_LINE << endl;
	cout << "> Checking cross section normalization for target '" << settings.targetNames[target_number] << "':" << endl;
	cout << "\tAnalytical integral: " << cs_int_analytical << " eV fm^2" << endl;
	cout << "\tIntegral of numerically calculated cross section: " << cs_int_numerical << " eV fm^2" << endl;
	cout << "\tDeviation: " << (cs_int_analytical - cs_int_numerical)/cs_int_analytical*100. << " %" << endl;
}
