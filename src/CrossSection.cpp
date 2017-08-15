#include "CrossSection.h"

#include <iostream>
#include <sstream>
#include <algorithm>

#include "omp.h"

#include "TGraph.h"
#include "TAxis.h"

using std::cout;
using std::endl;
using std::upper_bound;
using std::max_element;

void CrossSection::breit_wigner(double (&energy_bins)[NBINS], vector<double> (&crosssection_bins), double e0, double gamma0, double gamma, double jj, double j0){

	double cs_max = PI*0.5*HBARC2/(e0*e0)*(2.*jj + 1.)/(2. * j0 + 1.)*gamma0*gamma;
	
	for(unsigned int j = 0; j < NBINS; ++j){
		crosssection_bins[NBINS + j] += cs_max / ((energy_bins[j] - e0)*(energy_bins[j] - e0) + 0.25*gamma*gamma);
	}
}

void CrossSection::maxwell_boltzmann(double (&energy_bins)[NBINS], vector<double> &velocity_bins, vector<double> &vdist_bins, vector<double> &params, double mass, double e0){

// Implementation with non-equidistant velocity bins. Here, velocity_bins[i] is the velocity that is needed to shift energy_bins[j] to energy_bins[i].
	double c1 = sqrt(mass*AtomicMassUnit/(2*PI*kB*params[0]));
      	double c2 = -1./pow(delta(params[0], mass), 2);
	double energy_ratio = 1.;

	// Find bin of resonance energy
	long int resonance_position = upper_bound(energy_bins, energy_bins + NBINS, e0) - energy_bins;

      	for(unsigned int i = 0; i < NBINS; ++i){
		energy_ratio = energy_bins[i]/energy_bins[resonance_position];
		velocity_bins[i] = (-2. + 2.*energy_ratio*energy_ratio)/(2. + 2.*energy_ratio*energy_ratio);
              	vdist_bins[i] = c1*exp(c2*velocity_bins[i]*velocity_bins[i]);
      	}
}

void CrossSection::maxwell_boltzmann_debye(double (&energy_bins)[NBINS], vector<double> &velocity_bins, vector<double> &vdist_bins, vector<double> &params, double mass, double e0){
	;
}

void CrossSection::dopplershift(double (&dopplercs_bins)[NBINS], double (&energy_bins)[NBINS], vector<vector<double> > &crosssection_bins, vector< vector<double> > &velocity_bins, vector<vector<double> > &vdist_bins, vector<double> &vdist_norm){

	for(unsigned int i = 0; i < crosssection_bins.size(); ++i){
		unsigned long int resonance_position = (unsigned long int) (max_element(crosssection_bins[i].begin(), crosssection_bins[i].end()) - crosssection_bins[i].begin() + 1);

		#pragma omp parallel for
		for(unsigned int j = 0; j < NBINS; ++j){
			for(unsigned int k = 0; k < NBINS - 1; ++k){
				dopplercs_bins[j] += vdist_norm[i]*vdist_bins[i][k]*crosssection_bins[i][j + resonance_position - k]*(velocity_bins[i][k + 1] - velocity_bins[i][k]);
			}
		}
	}
}

void CrossSection::maxwell_boltzmann_approximation(double (&dopplercs_bins)[NBINS], double (&energy_bins)[NBINS], vector< vector<double> > &velocity_bins, vector< vector<double> > &vdist_bins, vector<double> &e0_list, vector<double> &gamma0_list, vector<double> &gamma_list, vector<double> &jj_list, double j0, vector<double> &params, double mass){

// Calculate velocity distribution as in maxwell_boltzmann

	for(unsigned int i = 0; i < e0_list.size(); ++i){
		velocity_bins.push_back(vector<double> (NBINS));
		vdist_bins.push_back(vector<double> (NBINS));
		double c1 = sqrt(mass*AtomicMassUnit/(2*PI*kB*params[0]));
		double c2 = -1./pow(delta(params[0], mass), 2);
		double energy_ratio = 1.;

		// Find bin of resonance energy
		long int resonance_position = upper_bound(energy_bins, energy_bins + NBINS, e0_list[i]) - energy_bins;

		for(unsigned int j = 0; j < NBINS; ++j){
			energy_ratio = energy_bins[j]/energy_bins[resonance_position];
			velocity_bins[i][j] = (-2. + 2.*energy_ratio*energy_ratio)/(2. + 2.*energy_ratio*energy_ratio);
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
			graphs[i]->SetPoint((Int_t) j, energy_bins[j], crosssection_bins[i][NBINS + j]);
			total_crosssection_bins[j] += crosssection_bins[i][NBINS + j];
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
			total_crosssection_bins[j] += crosssection_bins[i][NBINS + j];
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
