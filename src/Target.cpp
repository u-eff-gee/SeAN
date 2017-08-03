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
};

void Target::calculateVelocityDistribution(double (&vdist_bins)[NBINS_V], string vdist_ID){
	if(vdist_ID == "absolute_zero"){
		for(int i = 0; i < NBINS_E; ++i){
			dopplercs_bins[i] = crosssection_bins[i];
		}
	}
	if(vdist_ID == "maxwell_boltzmann"){
		crossSection->maxwell_boltzmann(vdist_bins, crosssection_bins, vDistParams);
	}
	if(vdist_ID == "maxwell_boltzmann_debye"){
		crossSection->maxwell_boltzmann_debye(energy_bins, crosssection_bins, vDistParams);
	}
}

void Target::plotCrossSection(double (&energy_bins)[NBINS_E]){

	stringstream filename;
	filename << target_name << "_cross_section.pdf";
	TCanvas *canvas = new TCanvas("canvas", target_name.c_str(), 0, 0, 800, 500);
	TLegend *legend = new TLegend(CS_PLOT_LEGEND_X1, CS_PLOT_LEGEND_Y1, CS_PLOT_LEGEND_X2, CS_PLOT_LEGEND_Y2);

	crossSection->plot(energy_bins, crosssection_bins, target_name, canvas, legend, "Cross section", false);

	legend->Draw();
	canvas->SaveAs(filename.str().c_str());
};

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
