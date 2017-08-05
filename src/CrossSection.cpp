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

void CrossSection::breit_wigner(double (&energy_bins)[NBINS], double (&crosssection)[NBINS], vector<double> &e0, vector<double> &gamma0, vector<double> &gamma, vector<double> &jj, double j0){

	double cs_max = 0.;

	for(unsigned int i = 0; i < e0.size(); ++i){
		cs_max = PI*0.5*HBARC2/(e0[i]*e0[i])*(2.*jj[i] + 1.)/(2. * j0 + 1.)*gamma0[i]*gamma[i];
		for(int j = 0; j < NBINS; ++j){
			crosssection[j] += cs_max / ((energy_bins[j] - e0[i])*(energy_bins[j] - e0[i]) + 0.25*gamma[i]*gamma[i]);
		}
	}
}

void CrossSection::maxwell_boltzmann(double (&energy_bins)[NBINS], vector<double> &velocity_bins, vector<double> &vdist_bins, vector<double> &params, double mass, double e0){
// Use the fact that the velocity distribution is symmetric and store only one half of the distribution

// Implementation with equidistant velocity bins.
//	double v_max = 3.*delta(params[0], mass);
//	double delta_v = v_max/NBINS;
//
//	double c1 = sqrt(mass*AtomicMassUnit/(2*PI*kB*params[0]));
//	double c2 = -1./pow(delta(params[0], mass), 2);
//
//
//	for(int i = 0; i < NBINS; ++i){
//		velocity_bins[i] = i*delta_v;
//		vdist_bins[i] = c1*exp(c2*velocity_bins[i]*velocity_bins[i]);
//	}

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

void CrossSection::dopplershift(double (&dopplercs_bins)[NBINS], double (&energy_bins)[NBINS], double (&crosssection_bins)[NBINS], vector<double> &velocity_bins, vector<double> &vdist_bins){
	for(int i = 0; i < NBINS; ++i){
		for(int j = 0; j < NBINS; ++j){
			dopplercs_bins[i] += vdist_bins[j]*crosssection_bins[j];
		}
	}
}

void CrossSection::plot_crosssection(double (&energy_bins)[NBINS], double (&crosssection_bins)[NBINS], string title, TCanvas* canvas, TLegend* legend, string legend_entry, bool add_to_existing_canvas){

	TGraph *graph = new TGraph(NBINS, energy_bins, crosssection_bins);
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

void CrossSection::plot_vdist(vector<double> &velocity_bins, vector<double> &vdist_bins, string title, TCanvas* canvas, TLegend* legend, string legend_entry, bool add_to_existing_canvas){

	TGraph *graph = new TGraph(NBINS);
	for(int i = 0; i < NBINS; ++i){
		graph->SetPoint(i, velocity_bins[i], vdist_bins[i]);
	}
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
