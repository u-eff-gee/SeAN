#include "Absorption.h"

#include <sstream>
#include <fstream>
#include <iostream>

using std::stringstream;
using std::ifstream;
using std::cout;
using std::endl;

void Absorption::read_massattenuation_NIST(double (&energy_bins)[NBINS], double (&massattenuation_bins)[NBINS], string massAttenuation_ID, double mass){
	
	stringstream filename;
	filename << MU_DIR << massAttenuation_ID << ".dat";

	string line;
	nbins_matt = 0;
	ifstream ifile;
	size_t start = 0;
	size_t stop = 0;

	ifile.open(filename.str().c_str());	

        if(!ifile.is_open()){
                cout << "Error: Absorption.cc: read_massattenuation(): File '" << filename.str() << "' not found." << endl;
		abort();
	}
        cout << "> Reading input file '" << filename.str() << "'" << endl;

	while(getline(ifile, line)){
		if(line.substr(0,1) == COMMENT)
			continue;
		
		start = 0;

		stop = line.find(NIST_SEPARATOR);

		matt[0].push_back(atof(line.substr(start, stop).c_str()));
		start = stop + NIST_SEPARATOR.length();
		stop = line.substr(start, line.length()).find(NIST_SEPARATOR);
		matt[1].push_back(atof(line.substr(start, stop).c_str()));

		++nbins_matt;
	}

	inter = new ROOT::Math::Interpolator(nbins_matt, ROOT::Math::Interpolation::kCSPLINE);
	inter->SetData(nbins_matt, &matt[0][0], &matt[1][0]);

	// Conversion from cm2/g to fm2/atom
	double conversion_factor = mass*AtomicMassUnitG*1.e26;

	for(int i = 0; i < NBINS; ++i){
		massattenuation_bins[i] = conversion_factor*inter->Eval(energy_bins[i]*0.000001);
	}
}

void Absorption::plot_massattenuation(double (&energy_bins)[NBINS], double (&massattenuation_bins)[NBINS], string title, TCanvas *canvas, TLegend* legend, string legend_entry){

	TGraph *graph = new TGraph(NBINS, energy_bins, massattenuation_bins);
	graph->SetName(title.c_str());
	graph->SetTitle(title.c_str());
	graph->GetXaxis()->SetTitle("Energy / eV");
	graph->GetYaxis()->SetTitle("#mu / fm^2/atom");
	graph->Draw();

	legend->AddEntry(graph->GetName(), legend_entry.c_str(), "l");
}

void Absorption::plot_total_massattenuation(string title, TCanvas *canvas, TLegend* legend, string legend_entry){

	TGraph *graph = new TGraph(nbins_matt, &matt[0][0], &matt[1][0]);
	
	graph->SetName(title.c_str());
	graph->SetTitle(title.c_str());
	graph->GetXaxis()->SetTitle("Energy / MeV");
	graph->GetYaxis()->SetTitle("#mu / cm^2/g");
	graph->Draw();

	legend->AddEntry(graph->GetName(), legend_entry.c_str(), "l");

}

