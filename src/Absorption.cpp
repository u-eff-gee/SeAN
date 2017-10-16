#include "Absorption.h"

#include <sstream>
#include <fstream>
#include <iostream>
#include <regex>

#include "TGraph2D.h"
#include "TH2D.h"
#include "TStyle.h"

using std::stringstream;
using std::ifstream;
using std::cout;
using std::endl;
using std::regex;
using std::regex_replace;

void Absorption::read_massattenuation_NIST(double (&energy_bins)[NBINS], double (&massattenuation_bins)[NBINS], string massAttenuation_ID, double mass){
	
	stringstream filename;
	filename << MU_DIR << massAttenuation_ID << ".dat";

	string line;
	nbins_matt = 0;
	ifstream ifile;

	ifile.open(filename.str().c_str());	

        if(!ifile.is_open()){
                cout << "Error: Absorption.cc: read_massattenuation(): File '" << filename.str() << "' not found." << endl;
		abort();
	}
        cout << "> Reading input file '" << filename.str() << "'" << endl;

	while(getline(ifile, line)){
		if(line.substr(0,1) == COMMENT)
			continue;

		// Ignore lines with x-ray resonances since they have the same energy value as the previous bins. Those steps can not be interpolated.
		if(regex_replace(line.substr(NIST_XRAY, NIST_XRAY_LENGTH), regex("\\s+"), "") != "")
			continue;	
		
		matt[0].push_back(atof(line.substr(NIST_ENERGY, NIST_ENERGY_LENGTH).c_str()));
		matt[1].push_back(atof(line.substr(NIST_MU, NIST_MU_LENGTH).c_str()));

		++nbins_matt;
	}

	inter = new ROOT::Math::Interpolator(nbins_matt - 1, ROOT::Math::Interpolation::kCSPLINE);
	inter->SetData(nbins_matt - 1, &matt[0][0], &matt[1][0]);

	// Conversion from cm2/g to fm2/atom
	double conversion_factor = mass*AtomicMassUnitG*1.e26;

	for(int i = 0; i < NBINS; ++i){
		massattenuation_bins[i] = conversion_factor*inter->Eval(energy_bins[i]*0.000001);
	}
}

void Absorption::const_beam(double (&energy_bins)[NBINS], double (&incident_beam_bins)[NBINS], vector<double> beamParams){

	for(int i = 0; i < NBINS; ++i)
		incident_beam_bins[i] = beamParams[0];	

}

void Absorption::gauss_beam(double (&energy_bins)[NBINS], double (&incident_beam_bins)[NBINS], vector<double> beamParams){

	double c1 = -1./(2.*beamParams[1]*beamParams[1]);

	for(int i = 0; i < NBINS; ++i)
		incident_beam_bins[i] = beamParams[2]*exp((energy_bins[i] - beamParams[0])*(energy_bins[i] - beamParams[0])*c1);

}

void Absorption::photon_flux_density(double (&dopplercs_bins)[NBINS], double (&massattenuation_bins)[NBINS], double (&z_bins)[NBINS_Z], double (&incident_beam_bins)[NBINS], double (&photon_flux_density_bins)[NBINS][NBINS_Z]){

	//#pragma omp parallel for
	for(int i = 0; i < NBINS; ++i){
		for(int j = 0; j < NBINS_Z; ++j){
			photon_flux_density_bins[i][j] = incident_beam_bins[i]*exp(-(dopplercs_bins[i] + massattenuation_bins[i])*z_bins[j]);
		}
	}
}

void Absorption::resonance_absorption_density(double (&dopplercs_bins)[NBINS], double (&photon_flux_density_bins)[NBINS][NBINS_Z], double (&resonance_absorption_density_bins)[NBINS][NBINS_Z]){

	//#pragma omp parallel for
	for(int i = 0; i < NBINS; ++i){
		for(int j = 0; j < NBINS_Z; ++j){
			resonance_absorption_density_bins[i][j] = dopplercs_bins[i]*photon_flux_density_bins[i][j];
		}
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

	TGraph *graph = new TGraph((Int_t) nbins_matt - 1, &matt[0][0], &matt[1][0]);
	
	graph->SetName(title.c_str());
	graph->SetTitle(title.c_str());
	graph->GetXaxis()->SetTitle("Energy / MeV");
	graph->GetYaxis()->SetTitle("#mu / cm^2/g");
	graph->Draw();

	legend->AddEntry(graph->GetName(), legend_entry.c_str(), "l");

}

void Absorption::plot_photon_flux_density(double (&energy_bins)[NBINS], double (&z_bins)[NBINS_Z], double (&photon_flux_density_bins)[NBINS][NBINS_Z], string title, TCanvas *canvas, TLegend* legend, string legend_entry){
	
	// Saveguard for the number of bins used for plotting phi. TGraph2D uses up a lot of memory if the dimension of the matrix for phi is too large.
	int nbins_plot = NBINS;
	int nbins_skip = 1;
	int nbins_z_plot = NBINS_Z;
	int nbins_z_skip = 1;

	if(NBINS > NBINS_PHI_MAX){
		nbins_skip = NBINS/NBINS_PHI_MAX + 1;
		nbins_plot = NBINS/nbins_skip;

		cout << "> Warning: Absorption.cpp: plot_photon_flux_density(): The number of energy bins (" << NBINS << ") is larger than the maximum number of bins that can be used for plotting the photon flux density (" << NBINS_PHI_MAX << "). Using only 1 out of " << nbins_skip << " bins." << endl;
		cout << "  If you are not satisfied with the quality of the plot, increase NBINS_PHI_MAX in Config.h and recompile the program." << endl;
	}
	if(NBINS > NBINS_Z_PHI_MAX){
		nbins_z_skip = NBINS_Z/NBINS_Z_PHI_MAX + 1;
		nbins_z_plot = NBINS_Z/nbins_z_skip;

		cout << "> Warning: Absorption.cpp: plot_photon_flux_density(): The number of z bins (" << NBINS_Z << ") is larger than the maximum number of bins that can be used for plotting the photon flux density (" << NBINS_Z_PHI_MAX << "). Using only 1 out of " << nbins_z_skip << " bins." << endl;
		cout << "  If you are not satisfied with the quality of the plot, increase NBINS_Z_PHI_MAX in Config.h and recompile the program." << endl;
	}

	TGraph2D *graph = new TGraph2D(nbins_plot*nbins_z_plot);	
	int npoint = 0;
	for(int i = 0; i < NBINS; i += nbins_skip){
		for(int j = 0; j < NBINS_Z; j += nbins_z_skip){
			graph->SetPoint(npoint, energy_bins[i], z_bins[j], photon_flux_density_bins[i][j]);
			++npoint;
		}
	}

	graph->SetName(title.c_str());
	graph->SetTitle(title.c_str());
	graph->GetHistogram()->GetXaxis()->SetTitle("Energy / eV");
	graph->GetHistogram()->GetXaxis()->SetTitleOffset(2.);
	graph->GetHistogram()->GetYaxis()->SetTitle("z / atoms/fm^{2}");
	graph->GetHistogram()->GetYaxis()->SetTitleOffset(2.);
	graph->GetHistogram()->GetZaxis()->SetTitle("#Phi (E, z)");
	gStyle->SetPalette(55);
	graph->Draw("surf1");

	canvas->SetTheta(45.);
	canvas->SetPhi(200.);

	legend->AddEntry(graph->GetName(), legend_entry.c_str(), "l");
}

void Absorption::plot_test_integration(double (&energy_bins)[NBINS], double (&z_bins)[NBINS_Z], double (&photon_flux_density_bins)[NBINS][NBINS_Z], string title, TCanvas *canvas, TLegend* legend, string legend_entry){
	
	// Saveguard for the number of bins used for plotting phi. TGraph2D uses up a lot of memory if the dimension of the matrix for phi is too large.
	int nbins_plot = NBINS;
	int nbins_skip = 1;
	int nbins_z_plot = NBINS_Z;
	int nbins_z_skip = 1;

	if(NBINS > NBINS_PHI_MAX){
		nbins_skip = NBINS/NBINS_PHI_MAX + 1;
		nbins_plot = NBINS/nbins_skip;

		cout << "> Warning: Absorption.cpp: plot_test_integration(): The number of x bins (" << NBINS << ") is larger than the maximum number of bins that can be used for plotting the integration test function f (" << NBINS_PHI_MAX << "). Using only 1 out of " << nbins_skip << " bins." << endl;
		cout << "  If you are not satisfied with the quality of the plot, increase NBINS_PHI_MAX in Config.h and recompile the program." << endl;
	}
	if(NBINS > NBINS_Z_PHI_MAX){
		nbins_z_skip = NBINS_Z/NBINS_Z_PHI_MAX + 1;
		nbins_z_plot = NBINS_Z/nbins_z_skip;

		cout << "> Warning: Absorption.cpp: plot_test_integration(): The number of y bins (" << NBINS_Z << ") is larger than the maximum number of bins that can be used for plotting the photon flux density (" << NBINS_Z_PHI_MAX << "). Using only 1 out of " << nbins_z_skip << " bins." << endl;
		cout << "  If you are not satisfied with the quality of the plot, increase NBINS_Z_PHI_MAX in Config.h and recompile the program." << endl;
	}

	TGraph2D *graph = new TGraph2D(nbins_plot*nbins_z_plot);	
	int npoint = 0;
	for(int i = 0; i < NBINS; i += nbins_skip){
		for(int j = 0; j < NBINS_Z; j += nbins_z_skip){
			graph->SetPoint(npoint, energy_bins[i], z_bins[j], photon_flux_density_bins[i][j]);
			++npoint;
		}
	}

	graph->SetName(title.c_str());
	graph->SetTitle(title.c_str());
	graph->GetHistogram()->GetXaxis()->SetTitle("x");
	graph->GetHistogram()->GetXaxis()->SetTitleOffset(2.);
	graph->GetHistogram()->GetYaxis()->SetTitle("y");
	graph->GetHistogram()->GetYaxis()->SetTitleOffset(2.);
	graph->GetHistogram()->GetZaxis()->SetTitle("z");
	gStyle->SetPalette(55);
	graph->Draw("surf1");

	canvas->SetTheta(45.);
	canvas->SetPhi(200.);

	legend->AddEntry(graph->GetName(), legend_entry.c_str(), "l");
}

void Absorption::plot_resonance_absorption_density(double (&energy_bins)[NBINS], double (&z_bins)[NBINS_Z], double (&resonance_absorption_density_bins)[NBINS][NBINS_Z], string title, TCanvas *canvas, TLegend* legend, string legend_entry){
	
	// Saveguard for the number of bins used for plotting phi. TGraph2D uses up a lot of memory if the dimension of the matrix for phi is too large.
	int nbins_plot = NBINS;
	int nbins_skip = 1;
	int nbins_z_plot = NBINS_Z;
	int nbins_z_skip = 1;

	if(NBINS > NBINS_ALPHA_MAX){
		nbins_skip = NBINS/NBINS_ALPHA_MAX + 1;
		nbins_plot = NBINS/nbins_skip;

		cout << "> Warning: Absorption.cpp: plot_resonance_absorption_density(): The number of energy bins (" << NBINS << ") is larger than the maximum number of bins that can be used for plotting the resonance absorption density (" << NBINS_ALPHA_MAX << "). Using only 1 out of " << nbins_skip << " bins." << endl;
		cout << "  If you are not satisfied with the quality of the plot, increase NBINS_ALPHA_MAX in Config.h and recompile the program." << endl;
	}
	if(NBINS > NBINS_Z_ALPHA_MAX){
		nbins_z_skip = NBINS_Z/NBINS_Z_ALPHA_MAX + 1;
		nbins_z_plot = NBINS_Z/nbins_z_skip;

		cout << "> Warning: Absorption.cpp: plot_resonance_absorption_density(): The number of z bins (" << NBINS_Z << ") is larger than the maximum number of bins that can be used for plotting the resonance absorption density (" << NBINS_Z_ALPHA_MAX << "). Using only 1 out of " << nbins_z_skip << " bins." << endl;
		cout << "  If you are not satisfied with the quality of the plot, increase NBINS_Z_ALPHA_MAX in Config.h and recompile the program." << endl;
	}

	TGraph2D *graph = new TGraph2D(nbins_plot*nbins_z_plot);	

	int npoint = 0;
	for(int i = 0; i < NBINS; i += nbins_skip){
		for(int j = 0; j < NBINS_Z; j += nbins_z_skip){
			graph->SetPoint(npoint, energy_bins[i], z_bins[j], resonance_absorption_density_bins[i][j]);
			++npoint;
		}
	}

	graph->SetName(title.c_str());
	graph->SetTitle(title.c_str());
	graph->GetHistogram()->GetXaxis()->SetTitle("Energy / eV");
	graph->GetHistogram()->GetXaxis()->SetTitleOffset(2.);
	graph->GetHistogram()->GetYaxis()->SetTitle("z / atoms/fm^{2}");
	graph->GetHistogram()->GetYaxis()->SetTitleOffset(2.);
	graph->GetHistogram()->GetZaxis()->SetTitle("#alpha (E, z) / fm^{2}");
	gStyle->SetPalette(77);
	graph->Draw("surf1");

	canvas->SetTheta(45.);
	canvas->SetPhi(200.);

	legend->AddEntry(graph->GetName(), legend_entry.c_str(), "l");
}
