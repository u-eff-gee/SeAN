#include "CrossSection.h"

#include <iostream>
#include <sstream>
#include <algorithm>

#include "TGraph.h"
#include "TAxis.h"

using std::cout;
using std::endl;
using std::stringstream;
using std::upper_bound;
using std::max_element;

void CrossSection::breit_wigner(double (&energy_bins)[NBINS], vector<double> (&crosssection_bins), double e0, double gamma0, double gamma, double jj, double j0){

	double cs_max = PI*0.5*HBARC2/(e0*e0)*(2.*jj + 1.)/(2. * j0 + 1.)*gamma0*gamma;
	
	for(int j = 0; j < NBINS; ++j){
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

      	for(int i = 0; i < NBINS; ++i){
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
		long int resonance_position = max_element(crosssection_bins[i].begin(), crosssection_bins[i].end()) - crosssection_bins[i].begin() + 1;

		for(int j = 0; j < NBINS; ++j){
			for(int k = 0; k < NBINS - 1; ++k){
				dopplercs_bins[j] += vdist_norm[i]*vdist_bins[i][k]*crosssection_bins[i][j + resonance_position - k]*(velocity_bins[i][k + 1] - velocity_bins[i][k]);
			}
		}
	}
}

void CrossSection::plot_crosssection(double (&energy_bins)[NBINS], vector< vector<double> > (&crosssection_bins), string title, TCanvas* canvas, TLegend* legend, string legend_entry){

	double total_crosssection_bins[NBINS] = {0.};
	vector< TGraph* > graphs;

	for(unsigned int i = 0; i < crosssection_bins.size(); ++i){
		graphs.push_back(new TGraph(NBINS));
		for(int j = 0; j < NBINS; ++j){
			graphs[i]->SetPoint(j, energy_bins[j], crosssection_bins[i][NBINS + j]);
			total_crosssection_bins[j] += crosssection_bins[i][NBINS + j];
		}

		graphs[i]->SetLineStyle(2);
	}

	TGraph* graph = new TGraph(NBINS, energy_bins, total_crosssection_bins);
	
	graph->SetName(title.c_str());
	graph->SetTitle(title.c_str());
	graph->GetXaxis()->SetTitle("Energy / eV");
	graph->GetYaxis()->SetTitle("#sigma / fm^2");
		
	legend->AddEntry(graph->GetName(), legend_entry.c_str(), "l");
	graph->Draw();

	for(unsigned int i = 0; i < graphs.size(); ++i)
		graphs[i]->Draw("same");
}

void CrossSection::plot_vdist(vector<double> &velocity_bins, vector<double> &vdist_bins, string title, TCanvas* canvas, TLegend* legend, string legend_entry){

	TGraph *graph = new TGraph(NBINS);
	for(int i = 0; i < NBINS; ++i){
		graph->SetPoint(i, velocity_bins[i], vdist_bins[i]);
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
		for(int j = 0; j < NBINS; ++j){
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
