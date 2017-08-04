#include "Target.h"

#include <iostream>
#include <sstream>

#include "TCanvas.h"
#include "TLegend.h"

using std::cout;
using std::endl;
using std::stringstream;

void Target::calculateCrossSection(double (&energy_bins)[NBINS_E]){
	crossSection->breit_wigner(energy_bins, crosssection_bins, e0_list, gamma0_list, gamma_list, jj_list, j0);	
}

void Target::plotCrossSection(double (&energy_bins)[NBINS_E]){

	stringstream filename;
	filename << target_name << "_cross_section.pdf";
	stringstream canvasname;
	canvasname << target_name << "_canvas";
	TCanvas *canvas = new TCanvas(canvasname.str().c_str(), target_name.c_str(), 0, 0, 800, 500);
	TLegend *legend = new TLegend(CS_PLOT_LEGEND_X1, CS_PLOT_LEGEND_Y1, CS_PLOT_LEGEND_X2, CS_PLOT_LEGEND_Y2);

	crossSection->plot_crosssection(energy_bins, crosssection_bins, target_name, canvas, legend, "Cross section", false);

	legend->Draw();
	canvas->SaveAs(filename.str().c_str());
	delete canvas;
}

void Target::calculateVelocityDistribution(double (&velocity_bins)[NBINS_V]){
	if(vDist_ID == "absolute_zero"){
		for(int i = 0; i < NBINS_E; ++i){
			dopplercs_bins[i] = crosssection_bins[i];
		}
	}
	if(vDist_ID == "maxwell_boltzmann"){
		crossSection->maxwell_boltzmann(velocity_bins, vdist_bins, vDistParams, mass);
	}
//	if(vDist_ID == "maxwell_boltzmann_debye"){
//		crossSection->maxwell_boltzmann_debye(velocity_bins, vdist_bins, vDistParams);
//	}
}

void Target::plotVelocityDistribution(double (&velocity_bins)[NBINS_V]){

	stringstream filename;
	filename << target_name << "_velocity_distribution.pdf";
	stringstream canvasname;
	canvasname << target_name << "_canvas";
	TCanvas *canvas = new TCanvas(canvasname.str().c_str(), target_name.c_str(), 0, 0, 800, 500);
	TLegend *legend = new TLegend(V_PLOT_LEGEND_X1, V_PLOT_LEGEND_Y1, V_PLOT_LEGEND_X2, V_PLOT_LEGEND_Y2);

	crossSection->plot_vdist(velocity_bins, vdist_bins, target_name, canvas, legend, "Velocity Distribution", false);

	legend->Draw();
	canvas->SaveAs(filename.str().c_str());
	delete canvas;
}

void Target::calculateDopplerShift(double (&velocity_bins)[NBINS_V], double (&energy_bins)[NBINS_E]){
	crossSection->dopplershift(dopplercs_bins, velocity_bins, vdist_bins, energy_bins, crosssection_bins);
}

void Target::print(){
	cout << "TARGET #" << target_number << ": '" << target_name << "'" << endl;
	cout << "GROUND STATE J = " << j0 << endl;
	cout << "RESONANCES:" << endl;
	cout << "E0/eV\tGAMMA0\tGAMMA\tJ" << endl;
	
	for(unsigned int i = 0; i < e0_list.size(); ++i)
		cout << e0_list[i] << "\t" << gamma0_list[i] << "\t" << gamma_list[i] << "\t" << jj_list[i] << endl;

	cout << "VELOCITY DISTRIBUTION = " << vDist_ID;

	if( vDistParams.size() ){
		cout << " ( ";
		for(unsigned int i = 0; i < vDistParams.size(); ++i)
			cout << vDistParams[i] << " ";

		cout << ")" << endl;
	} else{
		cout << endl;
	}

	cout << "MASS = " << mass << " u" << endl;
	cout << "MASS ATTENUATION = " << massAttenuation_ID << endl;
	cout << "TARGET THICKNESS = " << z << " atoms / fm^2" << endl;
}
