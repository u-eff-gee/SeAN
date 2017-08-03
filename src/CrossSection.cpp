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

void CrossSection::plot(double (&energy_bins)[NBINS_E], double (&crosssection_bins)[NBINS_E], string title, TCanvas* canvas, TLegend* legend, string legend_entry, bool add_to_existing_canvas){

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
