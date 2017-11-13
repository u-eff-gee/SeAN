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

void CrossSection::integration_input(vector< vector<double> > &crosssection_bins, vector< vector<double> > &vdist_bins){
	for(unsigned int i = 0; i < crosssection_bins.size(); ++i){
		pconv_crosssection_bins.push_back(vector<double>(3*NBINS, 0.));

		for(unsigned int j = 0; j < NBINS; ++j){
			pconv_crosssection_bins[i][j + NBINS] = crosssection_bins[i][j];
		}
	}
}

void CrossSection::fft_input(vector<double> &energy_bins, vector< vector<double> > &crosssection_bins, vector< vector<double> > &vdist_bins, vector<double> e0_list){
	double de = energy_bins[1] - energy_bins[0];
	double energy_ratio_squared = 1.;

	for(unsigned int i = 0; i < crosssection_bins.size(); ++i){
		pconv_crosssection_bins.push_back(vector<double>(NBINS, 0.));	
		pconv_vdist_bins.push_back(vector<double>(NBINS, 0.));	

		for(unsigned int j = 0; j < NBINS; ++j){
			pconv_crosssection_bins[i][j] = crosssection_bins[i][j]*de;

			energy_ratio_squared = energy_bins[j]*energy_bins[j]/(e0_list[i]*e0_list[i]);	
			pconv_vdist_bins[i][j] = -vdist_bins[i][j]*(-4.*energy_ratio_squared/e0_list[i])/((1. + energy_ratio_squared)*(1. + energy_ratio_squared));
		}
	}
}

void CrossSection::dopplershift(vector<double> &dopplercs_bins, vector<double> &energy_bins, vector<vector<double> > &crosssection_bins, vector< vector<double> > &velocity_bins, vector<vector<double> > &vdist_bins, vector<double> &vdist_norm, vector<double> &e0_list){

	// Integrate using trapezoidal rule
	
	double resonance_energy_squared = 1.;
	
	for(unsigned int n = 0; n < pconv_crosssection_bins.size(); ++n){
		long unsigned int resonance_position = (unsigned long int) (upper_bound(energy_bins.begin(), energy_bins.end(), e0_list[n]) - energy_bins.begin());
		resonance_energy_squared = e0_list[n]*e0_list[n];

		for(unsigned int i = 0; i < NBINS; ++i){
			for(unsigned int j = 0; j < NBINS - 1; ++j){
				dopplercs_bins[i] += 0.5*(
						vdist_bins[n][j]*pconv_crosssection_bins[n][NBINS + i + resonance_position - j]*resonance_energy_squared/(energy_bins[j]*energy_bins[j]) 
						+ vdist_bins[n][j + 1]*pconv_crosssection_bins[n][NBINS + i + resonance_position - (j + 1)]*resonance_energy_squared/(energy_bins[j + 1]*energy_bins[j + 1]))
					*(velocity_bins[n][j] - velocity_bins[n][j + 1]); 
			}
		}
	}
}

void CrossSection::dopplershiftFFT(vector<double> &dopplercs_bins, vector<double> &energy_bins, vector<vector<double> > &crosssection_bins, vector< vector<double> > &velocity_bins, vector<vector<double> > &vdist_bins, vector<double> &vdist_norm, vector<unsigned int> &vdist_centroid){

	fftw_plan vdist_plan, crosssection_plan, product_fft_plan;
	fftw_complex vdist_fft[NBINS/2 + 1] = {{0.}};
	fftw_complex crosssection_fft[NBINS/2 + 1] = {{0.}};
	fftw_complex product_fft[NBINS/2 + 1] = {{0.}};
	double convolution[NBINS] = {0.};

	for(unsigned int i = 0; i < crosssection_bins.size(); ++i){

		vdist_plan = fftw_plan_dft_r2c_1d(NBINS, &pconv_vdist_bins[i][0], vdist_fft, FFTW_ESTIMATE);
		crosssection_plan = fftw_plan_dft_r2c_1d(NBINS, &pconv_crosssection_bins[i][0], crosssection_fft, FFTW_ESTIMATE);

		// Calculate Fourier transform of both arrays
		fftw_execute(vdist_plan);
		fftw_execute(crosssection_plan);

		// Multiply the (complex) Fourier - transformed arrays
		for(unsigned int j = 0; j < NBINS / 2 + 1; ++j){
			product_fft[j][0] = crosssection_fft[j][0]*vdist_fft[j][0] - crosssection_fft[j][1]*vdist_fft[j][1];
			product_fft[j][1] = crosssection_fft[j][0]*vdist_fft[j][1] + crosssection_fft[j][1]*vdist_fft[j][0];
		}

		// Transform the product back
		product_fft_plan = fftw_plan_dft_c2r_1d(NBINS, product_fft, convolution, FFTW_ESTIMATE);
		fftw_execute(product_fft_plan);

		// Add the i-th convoluted cross section to the total cross section, observing that the transformation back changed the order of the elements in the array and leaves them scaled by NBINS.
		unsigned int shift = vdist_centroid[i];
		double inverse_norm = (double) 1./(vdist_norm[i]*NBINS);
		shift = NBINS - shift;
		for(unsigned int j = 0; j < NBINS; ++j){
			if(j + shift < NBINS){
				dopplercs_bins[j + shift] += convolution[j]*inverse_norm;
			} else{
				dopplercs_bins[j - (NBINS - shift)] += convolution[j]*inverse_norm;
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

void CrossSection::plot_crosssection(vector<double> &energy_bins, vector< vector<double> > (&crosssection_bins), string title, TCanvas* canvas, TLegend* legend, string legend_entry){

	double total_crosssection_bins[NBINS] = {0.};
	vector< TGraph* > graphs;

	for(unsigned int i = 0; i < crosssection_bins.size(); ++i){
		graphs.push_back(new TGraph(NBINS));
		for(unsigned int j = 0; j < NBINS; ++j){
			graphs[i]->SetPoint((Int_t) j, energy_bins[j], crosssection_bins[i][j]);
			total_crosssection_bins[j] += crosssection_bins[i][j];
		}

		graphs[i]->SetLineStyle(2);
	}

	TGraph* graph = new TGraph(NBINS, &energy_bins[0], total_crosssection_bins);
	
	graph->SetName(title.c_str());
	graph->SetTitle(title.c_str());
	graph->GetXaxis()->SetTitle("Energy / eV");
	graph->GetYaxis()->SetTitle("#sigma / fm^{2}");
		
	legend->AddEntry(graph->GetName(), legend_entry.c_str(), "l");
	graph->Draw();

	for(unsigned int i = 0; i < graphs.size(); ++i)
		graphs[i]->Draw("same");
}

void CrossSection::plot_vdist(vector<double> &velocity_bins, vector<double> &vdist_bins, string title, TCanvas* canvas, TLegend* legend, string legend_entry){

	TGraph *graph = new TGraph(NBINS);
	for(unsigned int i = 0; i < NBINS; ++i){
		graph->SetPoint((Int_t) i, velocity_bins[i], vdist_bins[i]);
	}
	graph->SetName(title.c_str());
	graph->SetTitle(title.c_str());
	graph->GetXaxis()->SetTitle("Velocity / c");
	graph->GetYaxis()->SetTitle("Velocity distribution");
	graph->Draw();

	legend->AddEntry(graph->GetName(), legend_entry.c_str(), "l");

}

void CrossSection::plot_dopplershift(vector<double> &energy_bins, vector< vector<double> > &crosssection_bins, vector<double> &dopplercs_bins, string title, TCanvas* canvas, TLegend* legend, string legend_entry){

	double total_crosssection_bins[NBINS] = {0.};

	for(unsigned int i = 0; i < crosssection_bins.size(); ++i){
		for(unsigned int j = 0; j < NBINS; ++j){
			total_crosssection_bins[j] += crosssection_bins[i][j];
		}
	}

	TGraph* graph = new TGraph(NBINS, &energy_bins[0], total_crosssection_bins);
	graph->SetName(title.c_str());
	graph->SetTitle(title.c_str());
	graph->GetXaxis()->SetTitle("Energy / eV");
	graph->GetYaxis()->SetTitle("#sigma / fm^2");
	graph->SetLineStyle(2);
	graph->Draw();

	legend->AddEntry(graph->GetName(), "Cross section", "l");

	graph = new TGraph(NBINS, &energy_bins[0], &dopplercs_bins[0]);
	graph->Draw("same");

	legend->AddEntry(graph->GetName(), legend_entry.c_str(), "l");
}
