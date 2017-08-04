#include "CrossSection.h"

#include <iostream>
#include <sstream>

#include "TGraph.h"
#include "TAxis.h"

using std::cout;
using std::endl;
using std::stringstream;

void CrossSection::breit_wigner(double (&energy_bins)[NBINS_E], double (&crosssection)[NBINS_E], vector<double> &e0, vector<double> &gamma0, vector<double> &gamma, vector<double> &jj, double j0){

	double cs_max = 0.;

	for(unsigned int i = 0; i < e0.size(); ++i){
		cs_max = PI*0.5*HBARC2/(e0[i]*e0[i])*(2.*jj[i] + 1.)/(2. * j0 + 1.)*gamma0[i]*gamma[i];
		for(int j = 0; j < NBINS_E; ++j){
			crosssection[j] += cs_max / ((energy_bins[j] - e0[i])*(energy_bins[j] - e0[i]) + 0.25*gamma[i]*gamma[i]);
		}
	}
}

void CrossSection::maxwell_boltzmann(double (&energy_bins)[NBINS_E], double (&velocity_bins)[NBINS_V], double (&vdist_bins)[NBINS_V], vector<double> &params, double mass){
// Use the fact that the velocity distribution is symmetric and store only one half of the distribution

// Implementation with equidistant velocity bins.
//	double v_max = 3.*delta(params[0], mass);
//	double delta_v = v_max/NBINS_V;
//
//	double c1 = sqrt(mass*AtomicMassUnit/(2*PI*kB*params[0]));
//	double c2 = -1./pow(delta(params[0], mass), 2);
//
//
//	for(int i = 0; i < NBINS_V; ++i){
//		velocity_bins[i] = i*delta_v;
//		vdist_bins[i] = c1*exp(c2*velocity_bins[i]*velocity_bins[i]);
//	}

// Implementation where velocity_bins[i] is the velocity that is needed to shift the resonance energy E0 to E0 + delta_E, where delta_E is the distance of the energy bins
	double c1 = sqrt(mass*AtomicMassUnit/(2*PI*kB*params[0]));
      	double c2 = -1./pow(delta(params[0], mass), 2);
	double energy_ratio = 1.;

      	for(int i = 0; i < NBINS_V; ++i){
		energy_ratio = energy_bins[i]/energy_bins[0];
		velocity_bins[i] = (-2. + 2.*energy_ratio*energy_ratio)/(2. + 2.*energy_ratio*energy_ratio);
              	vdist_bins[i] = c1*exp(c2*velocity_bins[i]*velocity_bins[i]);
      	}
}

void CrossSection::maxwell_boltzmann_debye(double (&energy_bins)[NBINS_E], double (&velocity_bins)[NBINS_V], double (&vdist_bins)[NBINS_V], vector<double> &params, double mass){
	;
}

void CrossSection::dopplershift(double (&dopplercs_bins)[NBINS_E], double (&energy_bins)[NBINS_E], double (&crosssection_bins)[NBINS_V], double (&velocity_bins)[NBINS_V], double (&vdist_bins)[NBINS_E]){
	for(int i = 0; i < NBINS_E; ++i){
		for(int j = 0; j < NBINS_V; ++j){
			dopplercs_bins[i] += 0.;
		}
	}
}

void CrossSection::plot_crosssection(double (&energy_bins)[NBINS_E], double (&crosssection_bins)[NBINS_E], string title, TCanvas* canvas, TLegend* legend, string legend_entry, bool add_to_existing_canvas){

	TGraph *graph = new TGraph(NBINS_E, energy_bins, crosssection_bins);
	graph->SetName(title.c_str());
	graph->SetTitle(title.c_str());
	graph->GetXaxis()->SetTitle("Energy / eV");
	graph->GetYaxis()->SetTitle("#sigma / fm^2");
	if(add_to_existing_canvas){
		graph->Draw("same");
	} else{
		graph->Draw();
	}

	legend->AddEntry(graph->GetName(), legend_entry.c_str(), "l");

}

void CrossSection::plot_vdist(double (&velocity_bins)[NBINS_V], double (&vdist_bins)[NBINS_V], string title, TCanvas* canvas, TLegend* legend, string legend_entry, bool add_to_existing_canvas){

	TGraph *graph = new TGraph(NBINS_V, velocity_bins, vdist_bins);
	graph->SetName(title.c_str());
	graph->SetTitle(title.c_str());
	graph->GetXaxis()->SetTitle("Velocity / c");
	graph->GetYaxis()->SetTitle("Velocity distribution");
	if(add_to_existing_canvas){
		graph->Draw("same");
	} else{
		graph->Draw();
	}

	legend->AddEntry(graph->GetName(), legend_entry.c_str(), "l");

}
