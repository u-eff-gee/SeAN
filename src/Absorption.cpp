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

//void Absorption::read_massattenuation_NIST(vector<double> &energy_bins, vector<double> &massattenuation_bins, string massAttenuation_ID, double mass){
//	
//	stringstream filename;
//	filename << "mass_attenuation" << massAttenuation_ID << ".dat";
//
//	string line;
//	nbins_matt = 0;
//	ifstream ifile;
//
//	ifile.open(filename.str().c_str());	
//
//        if(!ifile.is_open()){
//                cout << "Error: Absorption.cc: read_massattenuation(): File '" << filename.str() << "' not found." << endl;
//		abort();
//	}
//        cout << "> Reading input file '" << filename.str() << "'" << endl;
//
//	while(getline(ifile, line)){
//		if(line.substr(0,1) == COMMENT)
//			continue;
//
//		// Ignore lines with x-ray resonances since they have the same energy value as the previous bins. Those steps can not be interpolated.
//		if(regex_replace(line.substr(NIST_XRAY, NIST_XRAY_LENGTH), regex("\\s+"), "") != "")
//			continue;	
//		
//		matt[0].push_back(atof(line.substr(NIST_ENERGY, NIST_ENERGY_LENGTH).c_str()));
//		matt[1].push_back(atof(line.substr(NIST_MU, NIST_MU_LENGTH).c_str()));
//
//		++nbins_matt;
//	}
//
//	inter = new ROOT::Math::Interpolator(nbins_matt - 1, ROOT::Math::Interpolation::kCSPLINE);
//	inter->SetData(nbins_matt - 1, &matt[0][0], &matt[1][0]);
//
//	// Conversion from cm2/g to fm2/atom
//	double conversion_factor = mass*AtomicMassUnitG*1.e26;
//
//	for(unsigned int i = 0; i < NBINS; ++i){
//		massattenuation_bins[i] = conversion_factor*inter->Eval(energy_bins[i]*0.000001);
//	}
//}

void Absorption::const_beam(const vector<double> &energy_bins, vector<double> &incident_beam_histogram){

	for(unsigned int i = 0; i < settings.nbins_e; ++i)
		incident_beam_histogram[i] = settings.incidentBeamParams[0];	

}

void Absorption::gauss_beam(const vector<double> &energy_bins, vector<double> &incident_beam_histogram){

	double c1 = -1./(2.*settings.incidentBeamParams[1]*settings.incidentBeamParams[1]);

	for(unsigned int i = 0; i < settings.nbins_e; ++i)
		incident_beam_histogram[i] = settings.incidentBeamParams[2]*exp((energy_bins[i] - settings.incidentBeamParams[0])*(energy_bins[i] - settings.incidentBeamParams[0])*c1);

}

void Absorption::arbitrary_beam(const vector<double> &energy_bins, vector<double> &incident_beam_histogram, const vector< vector<double> > &incident_beam_file){

	// Interpolate data from file	
	ROOT::Math::Interpolator inter((unsigned int) incident_beam_file[0].size(), ROOT::Math::Interpolation::kCSPLINE);
	inter.SetData((unsigned int) incident_beam_file[0].size(), &incident_beam_file[0][0], &incident_beam_file[1][0]);

	// Evaluate at the given velocity bins
	for(unsigned int i = 0; i < incident_beam_histogram.size(); ++i){
		incident_beam_histogram[i] = inter.Eval(energy_bins[i]);		
	}
}

//void Absorption::photon_flux_density(vector<double> &dopplercs_bins, vector<double> &massattenuation_bins, vector<double> &z_bins, vector<double> &incident_beam_bins, vector<vector<double> > &photon_flux_density_bins){
//
//	//#pragma omp parallel for
//	for(unsigned int i = 0; i < NBINS; ++i){
//		for(unsigned int j = 0; j < NBINS_Z; ++j){
//			photon_flux_density_bins[i][j] = incident_beam_bins[i]*exp(-(dopplercs_bins[i] + massattenuation_bins[i])*z_bins[j]);
//		}
//	}
//}
//
//void Absorption::resonance_absorption_density(vector<double> &dopplercs_bins, vector<vector<double> > &photon_flux_density_bins, vector<vector<double> > &resonance_absorption_density_bins){
//
//	//#pragma omp parallel for
//	for(unsigned int i = 0; i < NBINS; ++i){
//		for(unsigned int j = 0; j < NBINS_Z; ++j){
//			resonance_absorption_density_bins[i][j] = dopplercs_bins[i]*photon_flux_density_bins[i][j];
//		}
//	}
//}
//
//void Absorption::plot_massattenuation(vector<double> &energy_bins, vector<double> &massattenuation_bins, string title, TCanvas *canvas, TLegend* legend, string legend_entry){
//
//	TGraph *graph = new TGraph(NBINS, &energy_bins[0], &massattenuation_bins[0]);
//	graph->SetName(title.c_str());
//	graph->SetTitle(title.c_str());
//	graph->GetXaxis()->SetTitle("Energy / eV");
//	graph->GetYaxis()->SetTitle("#mu / fm^2/atom");
//	graph->Draw();
//
//	legend->AddEntry(graph->GetName(), legend_entry.c_str(), "l");
//}
//
//void Absorption::plot_total_massattenuation(string title, TCanvas *canvas, TLegend* legend, string legend_entry){
//
//	TGraph *graph = new TGraph((Int_t) nbins_matt - 1, &matt[0][0], &matt[1][0]);
//	
//	graph->SetName(title.c_str());
//	graph->SetTitle(title.c_str());
//	graph->GetXaxis()->SetTitle("Energy / MeV");
//	graph->GetYaxis()->SetTitle("#mu / cm^2/g");
//	graph->Draw();
//
//	legend->AddEntry(graph->GetName(), legend_entry.c_str(), "l");
//
//}
//
//void Absorption::plot_photon_flux_density(vector<double> &energy_bins, vector<double> &z_bins, vector<vector<double> > &photon_flux_density_bins, string title, TCanvas *canvas, TLegend* legend, string legend_entry){
//	
//	// Saveguard for the number of bins used for plotting phi. TGraph2D uses up a lot of memory if the dimension of the matrix for phi is too large.
//	unsigned int nbins_plot = NBINS;
//	unsigned int nbins_skip = 1;
//	unsigned int nbins_z_plot = NBINS_Z;
//	unsigned int nbins_z_skip = 1;
//
//	if(NBINS > NBINS_PHI_MAX){
//		nbins_skip = NBINS/NBINS_PHI_MAX + 1;
//		nbins_plot = NBINS/nbins_skip;
//
//		cout << "> Warning: Absorption.cpp: plot_photon_flux_density(): The number of energy bins (" << NBINS << ") is larger than the maximum number of bins that can be used for plotting the photon flux density (" << NBINS_PHI_MAX << "). Using only 1 out of " << nbins_skip << " bins." << endl;
//		cout << "  If you are not satisfied with the quality of the plot, increase NBINS_PHI_MAX in Config.h and recompile the program." << endl;
//	}
//	if(NBINS > NBINS_Z_PHI_MAX){
//		nbins_z_skip = NBINS_Z/NBINS_Z_PHI_MAX + 1;
//		nbins_z_plot = NBINS_Z/nbins_z_skip;
//
//		cout << "> Warning: Absorption.cpp: plot_photon_flux_density(): The number of z bins (" << NBINS_Z << ") is larger than the maximum number of bins that can be used for plotting the photon flux density (" << NBINS_Z_PHI_MAX << "). Using only 1 out of " << nbins_z_skip << " bins." << endl;
//		cout << "  If you are not satisfied with the quality of the plot, increase NBINS_Z_PHI_MAX in Config.h and recompile the program." << endl;
//	}
//
//	TGraph2D *graph = new TGraph2D((int) (nbins_plot*nbins_z_plot));	
//	int npoint = 0;
//	for(unsigned int i = 0; i < NBINS; i += nbins_skip){
//		for(unsigned int j = 0; j < NBINS_Z; j += nbins_z_skip){
//			graph->SetPoint(npoint, energy_bins[i], z_bins[j], photon_flux_density_bins[i][j]);
//			++npoint;
//		}
//	}
//
//	graph->SetName(title.c_str());
//	graph->SetTitle(title.c_str());
//	graph->GetHistogram()->GetXaxis()->SetTitle("Energy / eV");
//	graph->GetHistogram()->GetXaxis()->SetTitleOffset(2.);
//	graph->GetHistogram()->GetYaxis()->SetTitle("z / atoms/fm^{2}");
//	graph->GetHistogram()->GetYaxis()->SetTitleOffset(2.);
//	graph->GetHistogram()->GetZaxis()->SetTitle("#Phi (E, z)");
//	gStyle->SetPalette(55);
//	graph->Draw("surf1");
//
//	canvas->SetTheta(45.);
//	canvas->SetPhi(200.);
//
//	legend->AddEntry(graph->GetName(), legend_entry.c_str(), "l");
//}
//
//void Absorption::plot_resonance_absorption_density(vector<double> &energy_bins, vector<double> &z_bins, vector<vector<double> > &resonance_absorption_density_bins, string title, TCanvas *canvas, TLegend* legend, string legend_entry){
//	
//	// Saveguard for the number of bins used for plotting phi. TGraph2D uses up a lot of memory if the dimension of the matrix for phi is too large.
//	unsigned int nbins_plot = NBINS;
//	unsigned int nbins_skip = 1;
//	unsigned int nbins_z_plot = NBINS_Z;
//	unsigned int nbins_z_skip = 1;
//
//	if(NBINS > NBINS_ALPHA_MAX){
//		nbins_skip = NBINS/NBINS_ALPHA_MAX + 1;
//		nbins_plot = NBINS/nbins_skip;
//
//		cout << "> Warning: Absorption.cpp: plot_resonance_absorption_density(): The number of energy bins (" << NBINS << ") is larger than the maximum number of bins that can be used for plotting the resonance absorption density (" << NBINS_ALPHA_MAX << "). Using only 1 out of " << nbins_skip << " bins." << endl;
//		cout << "  If you are not satisfied with the quality of the plot, increase NBINS_ALPHA_MAX in Config.h and recompile the program." << endl;
//	}
//	if(NBINS > NBINS_Z_ALPHA_MAX){
//		nbins_z_skip = NBINS_Z/NBINS_Z_ALPHA_MAX + 1;
//		nbins_z_plot = NBINS_Z/nbins_z_skip;
//
//		cout << "> Warning: Absorption.cpp: plot_resonance_absorption_density(): The number of z bins (" << NBINS_Z << ") is larger than the maximum number of bins that can be used for plotting the resonance absorption density (" << NBINS_Z_ALPHA_MAX << "). Using only 1 out of " << nbins_z_skip << " bins." << endl;
//		cout << "  If you are not satisfied with the quality of the plot, increase NBINS_Z_ALPHA_MAX in Config.h and recompile the program." << endl;
//	}
//
//	TGraph2D *graph = new TGraph2D((int) (nbins_plot*nbins_z_plot));	
//
//	int npoint = 0;
//	for(unsigned int i = 0; i < NBINS; i += nbins_skip){
//		for(unsigned int j = 0; j < NBINS_Z; j += nbins_z_skip){
//			graph->SetPoint(npoint, energy_bins[i], z_bins[j], resonance_absorption_density_bins[i][j]);
//			++npoint;
//		}
//	}
//
//	graph->SetName(title.c_str());
//	graph->SetTitle(title.c_str());
//	graph->GetHistogram()->GetXaxis()->SetTitle("Energy / eV");
//	graph->GetHistogram()->GetXaxis()->SetTitleOffset(2.);
//	graph->GetHistogram()->GetYaxis()->SetTitle("z / atoms/fm^{2}");
//	graph->GetHistogram()->GetYaxis()->SetTitleOffset(2.);
//	graph->GetHistogram()->GetZaxis()->SetTitle("#alpha (E, z) / fm^{2}");
//	gStyle->SetPalette(77);
//	graph->Draw("surf1");
//
//	canvas->SetTheta(45.);
//	canvas->SetPhi(200.);
//
//	legend->AddEntry(graph->GetName(), legend_entry.c_str(), "l");
//}
