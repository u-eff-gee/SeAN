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

void CrossSection::breit_wigner(double (&energy_bins)[NBINS], vector<double> (&crosssection_bins), double e0, double gamma0, double gamma, double jj, double j0){

	double cs_max = PI*0.5*HBARC2/(e0*e0)*(2.*jj + 1.)/(2.*j0 + 1.)*gamma0*gamma;
	
	for(unsigned int j = 0; j < NBINS; ++j){
		crosssection_bins[j] += cs_max / ((energy_bins[j] - e0)*(energy_bins[j] - e0) + 0.25*gamma*gamma);
	}
}

void CrossSection::calculateVelocityBins(double (&energy_bins)[NBINS], vector<double> &velocity_bins, double e0){

	// SeAN uses non-equidistant velocity bins. Here, velocity_bins[i] is the velocity that is needed to shift energy_bins[j] to energy_bins[i].
	double energy_ratio = 1.;

	// Find bin of resonance energy
//	long int resonance_position = upper_bound(energy_bins, energy_bins + NBINS, e0) - energy_bins;
//	double resonance_energy = energy_bins[resonance_position];

	for(unsigned int i = 0; i < NBINS; ++i){
		energy_ratio = energy_bins[i]/e0;
		velocity_bins[i] = (1. - energy_ratio*energy_ratio)/(1. + energy_ratio*energy_ratio);
	}
}

void CrossSection::maxwell_boltzmann(double (&energy_bins)[NBINS], vector<double> &velocity_bins, vector<double> &vdist_bins, vector<double> &params, double mass, double e0){

	double c1 = sqrt(mass*AtomicMassUnit/(2*PI*kB*params[0]));
      	double c2 = -1./pow(delta(params[0], mass), 2);

      	for(unsigned int i = 0; i < NBINS; ++i){
              	vdist_bins[i] = c1*exp(c2*velocity_bins[i]*velocity_bins[i]);
      	}
}

void CrossSection::maxwell_boltzmann_debye(double (&energy_bins)[NBINS], vector<double> &velocity_bins, vector<double> &vdist_bins, vector<double> &params, double mass, double e0){
	;
}

void CrossSection::fft_input(double (&energy_bins)[NBINS], vector< vector<double> > &crosssection_bins, vector< vector<double> > &vdist_bins, vector<double> e0_list){
	double de = energy_bins[1] - energy_bins[0];
	double energy_ratio_squared = 1.;

	for(unsigned int i = 0; i < crosssection_bins.size(); ++i){
		fft_crosssection_bins.push_back(vector<double>(NBINS, 0.));	
		fft_vdist_bins.push_back(vector<double>(NBINS, 0.));	

		for(unsigned int j = 0; j < NBINS; ++j){
			fft_crosssection_bins[i][j] = crosssection_bins[i][j]*de;

			energy_ratio_squared = energy_bins[j]*energy_bins[j]/(e0_list[i]*e0_list[i]);	
			fft_vdist_bins[i][j] = -vdist_bins[i][j]*(-4.*energy_ratio_squared/e0_list[i])/((1. + energy_ratio_squared)*(1. + energy_ratio_squared));
		}
	}
}

void CrossSection::dopplershift(double (&dopplercs_bins)[NBINS], double (&energy_bins)[NBINS], vector<vector<double> > &crosssection_bins, vector< vector<double> > &velocity_bins, vector<vector<double> > &vdist_bins, vector<double> &vdist_norm){

	fftw_plan vdist_plan, crosssection_plan, product_fft_plan;
	fftw_complex vdist_fft[NBINS/2 + 1] = {{0.}};
	fftw_complex crosssection_fft[NBINS/2 + 1] = {{0.}};
	fftw_complex product_fft[NBINS/2 + 1] = {{0.}};
	double convolution[NBINS] = {0.};

	for(unsigned int i = 0; i < crosssection_bins.size(); ++i){

		vdist_plan = fftw_plan_dft_r2c_1d(NBINS, &fft_vdist_bins[i][0], vdist_fft, FFTW_ESTIMATE);
		crosssection_plan = fftw_plan_dft_r2c_1d(NBINS, &fft_crosssection_bins[i][0], crosssection_fft, FFTW_ESTIMATE);

		// Calculate Fourier transform of both arrays
		fftw_execute(vdist_plan);
		fftw_execute(crosssection_plan);

		// Multiply the (complex) Fourier - transformed arrays
		for(unsigned int j = 0; j < NBINS / 2 + 1; ++j){
			product_fft[j][0] = crosssection_fft[j][0]*vdist_fft[j][0] - crosssection_fft[j][1]*vdist_fft[j][1];
			product_fft[j][1] = crosssection_fft[j][0]*vdist_fft[j][1] + crosssection_fft[j][1]*vdist_fft[j][0];

			//cout << crosssection_fft[j][0] << " + i * " << crosssection_fft[j][1] << endl;
			//cout << vdist_fft[j][0] << " + i * " << vdist_fft[j][1] << endl;
			//cout << product_fft[j][0] << " + i * " << product_fft[j][1] << endl;
		}

		// Transform the product back
		product_fft_plan = fftw_plan_dft_c2r_1d(NBINS, product_fft, convolution, FFTW_ESTIMATE);
		fftw_execute(product_fft_plan);

		// Add the i-th convoluted cross section to the total cross section, observing that the transformation back changed the order of the elements in the array and leaves them scaled by NBINS.
		for(unsigned int j = 0; j < NBINS; ++j){
			if(j < NBINS/2){
				dopplercs_bins[j] += (double) convolution[NBINS/2 - j - 1]/NBINS/fabs(vdist_norm[i]);
			} else{
				dopplercs_bins[j] += (double) convolution[NBINS - (j - NBINS/2) - 1]/NBINS/fabs(vdist_norm[i]);
			}
		}
	}
}

void CrossSection::maxwell_boltzmann_approximation(double (&dopplercs_bins)[NBINS], double (&energy_bins)[NBINS], vector< vector<double> > &velocity_bins, vector< vector<double> > &vdist_bins, vector<double> &e0_list, vector<double> &gamma0_list, vector<double> &gamma_list, vector<double> &jj_list, double j0, vector<double> &params, double mass){

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
		
		for(int j = 0; j < NBINS; ++j){
			dopplercs_bins[j] += cs_max*exp(-(energy_bins[j] - e0_list[i])*(energy_bins[j] - e0_list[i])/(doppler_width*doppler_width));
		}
	}
}

void CrossSection::plot_crosssection(double (&energy_bins)[NBINS], vector< vector<double> > (&crosssection_bins), string title, TCanvas* canvas, TLegend* legend, string legend_entry){

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

	TGraph* graph = new TGraph(NBINS, energy_bins, total_crosssection_bins);
	
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

void CrossSection::plot_dopplershift(double (&energy_bins)[NBINS], vector< vector<double> > &crosssection_bins, double (&dopplercs_bins)[NBINS], string title, TCanvas* canvas, TLegend* legend, string legend_entry){

	double total_crosssection_bins[NBINS] = {0.};

	for(unsigned int i = 0; i < crosssection_bins.size(); ++i){
		for(unsigned int j = 0; j < NBINS; ++j){
			total_crosssection_bins[j] += crosssection_bins[i][j];
		}
	}

	TGraph* graph = new TGraph(NBINS, energy_bins, total_crosssection_bins);
	graph->SetName(title.c_str());
	graph->SetTitle(title.c_str());
	graph->GetXaxis()->SetTitle("Energy / eV");
	graph->GetYaxis()->SetTitle("#sigma / fm^2");
	graph->SetLineStyle(2);
	graph->Draw();

	legend->AddEntry(graph->GetName(), "Cross section", "l");

	graph = new TGraph(NBINS, energy_bins, dopplercs_bins);
	graph->Draw("same");

	legend->AddEntry(graph->GetName(), legend_entry.c_str(), "l");
}
