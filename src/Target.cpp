#include "Target.h"

#include <iostream>
#include <sstream>

#include "TCanvas.h"
#include "TLegend.h"

using std::cout;
using std::endl;
using std::stringstream;

void Target::calculateCrossSection(double (&energy_bins)[NBINS]){
	for(unsigned int i = 0; i < e0_list.size(); ++i){
		// The vector of cross section bins is 3 times as large as the required energy range, because 
		crosssection_bins.push_back(vector <double>(3*NBINS, 0.));
		crossSection->breit_wigner(energy_bins, crosssection_bins[i], e0_list[i], gamma0_list[i], gamma_list[i], jj_list[i], j0);
	}
}

void Target::plotCrossSection(double (&energy_bins)[NBINS]){

	stringstream filename;
	filename << target_name << "_cross_section.pdf";
	stringstream canvasname;
	canvasname << target_name << "_canvas";
	TCanvas *canvas = new TCanvas(canvasname.str().c_str(), target_name.c_str(), 0, 0, 800, 500);
	TLegend *legend = new TLegend(CS_PLOT_LEGEND_X1, CS_PLOT_LEGEND_Y1, CS_PLOT_LEGEND_X2, CS_PLOT_LEGEND_Y2);

	crossSection->plot_crosssection(energy_bins, crosssection_bins, target_name, canvas, legend, "Cross section");

	legend->Draw();
	canvas->SaveAs(filename.str().c_str());
	delete canvas;
}

void Target::calculateVelocityDistribution(double (&energy_bins)[NBINS]){
	if(vDist_ID == "absolute_zero"){
		for(unsigned int i = 0; i < NBINS; ++i){
			for(unsigned int j = 0; j < crosssection_bins.size(); ++j){
				dopplercs_bins[i] = crosssection_bins[i][j];
			}
		}
	}

	if(vDist_ID == "maxwell_boltzmann"){
		for(unsigned int i = 0; i < e0_list.size(); ++i){
			velocity_bins.push_back(vector<double> (NBINS));
			vdist_bins.push_back(vector<double> (NBINS));
			crossSection->maxwell_boltzmann(energy_bins, velocity_bins[i], vdist_bins[i], vDistParams, mass, e0_list[i]);
			vdist_norm.push_back(normalizeVDist(i));
		}
	}

	if(vDist_ID == "maxwell_boltzmann_approximation"){
		crossSection->maxwell_boltzmann_approximation(dopplercs_bins, energy_bins, velocity_bins, vdist_bins, e0_list, gamma0_list, gamma0_list, jj_list, j0, vDistParams, mass);
	}

//	if(vDist_ID == "maxwell_boltzmann_debye"){
//		crossSection->maxwell_boltzmann_debye(velocity_bins, vdist_bins, vDistParams);
//	}
}

void Target::calculateIncidentBeam(double (&energy_bins)[NBINS], string beam_ID, vector<double> beamParams){
	if(beam_ID == "const")
		absorption->const_beam(energy_bins, incident_beam_bins, beamParams);

	if(beam_ID == "gauss")
		absorption->gauss_beam(energy_bins, incident_beam_bins, beamParams);
}

void Target::calculateZBins(){
	
	double delta_z = z/NBINS_Z;

	for(int i = 0; i < NBINS_Z; ++i){
		z_bins[i] = i*delta_z;
	}
}

void Target::calculateDopplerShift(double (&energy_bins)[NBINS]){
	if(vDist_ID != "maxwell_boltzmann_approximation"){
		crossSection->dopplershift(dopplercs_bins, energy_bins, crosssection_bins, velocity_bins, vdist_bins, vdist_norm);
	}
}

void Target::calculateMassAttenuation(double (&energy_bins)[NBINS]){
	if(massAttenuation_ID == "0"){
	;} else{
		absorption->read_massattenuation_NIST(energy_bins, massattenuation_bins, massAttenuation_ID, mass);
	}
}

void Target::setIncidentBeam(double &trans_beam_bins){
	for(int i = 0; i < NBINS; ++i)
		incident_beam_bins[i] = (&trans_beam_bins)[i];
}

void Target::calculatePhotonFluxDensity(){
	absorption->photon_flux_density(dopplercs_bins, massattenuation_bins, z_bins, incident_beam_bins, photon_flux_density_bins);
}

void Target::calculateTransmittedBeam(){
	for(int i = 0; i < NBINS; ++i)
		transmitted_beam_bins[i] = photon_flux_density_bins[i][NBINS_Z - 1];
}

void Target::calculateResonanceAbsorptionDensity(){
	absorption->resonance_absorption_density(dopplercs_bins, photon_flux_density_bins, resonance_absorption_density_bins);
}

double Target::integrateEZHistogram(double (&energy_bins)[NBINS], double (&z_bins)[NBINS_Z], double (&ezhist)[NBINS][NBINS_Z]){

	// Area of a bin in 2D plane
	double bin_area = (energy_bins[1] - energy_bins[0])*(z_bins[1] - z_bins[0]);

	// Crude implementation of the Riemann integral as a sum of bin contents times their dimension in energy- and z-direction. Since there are only NBINS-1 spaces between NBINS bins, leave out the last bin in each loop.
	double integral = 0.;

	for(int i = 0; i < NBINS - 1; ++i){
		for(int j = 0; j < NBINS_Z - 1; ++j){
			integral += bin_area*ezhist[i][j]; 
		}
	}

	// Implementation that starts at the 1st and not the 0th bin to check the validity of this approximation
//	double integral = 0.;
//	for(int i = 1; i < NBINS; ++i){
//		for(int j = 1; j < NBINS_Z; ++j){
//			integral += bin_area*ezhist[i][j]; 
//		}
//	}

	return integral;
}

void Target::testIntegration(double (&energy_bins)[NBINS], vector<double> beamParams){

	cout << "> Test of integration ..." << endl;
	cout << "[EMIN, EMAX] = [" << energy_bins[0] << ", " << energy_bins[NBINS - 1] << "]" << endl;
	cout << "c\tx0\ty0\tsigma" << endl;

	for(unsigned int i = 0; i < gamma0_list.size(); ++i){
		cout << beamParams[i] << "\t" << e0_list[2*i] << "\t" << e0_list[2*i + 1] << "\t" << gamma0_list[i] << endl;
	}

	double denominator = 1.;

	for(int i = 0; i < NBINS; ++i){
		for(int j = 0; j < NBINS_Z; ++j){
			photon_flux_density_bins[i][j] = 0.;
			// Use photon_flux_density_bins to store the test function
			for(unsigned int k = 0; k < beamParams.size(); ++k){
				denominator = 1./(2.*gamma0_list[k]*gamma0_list[k]);
				photon_flux_density_bins[i][j] += beamParams[k]*exp(-denominator*(energy_bins[i] - e0_list[2*k])*(energy_bins[i] - e0_list[2*k]))*exp(-denominator*(energy_bins[j] - e0_list[2*k + 1])*(energy_bins[j] - e0_list[2*k + 1]));
			}
		}
	}

	//plotTestIntegration(energy_bins);

	double result = integrateEZHistogram(energy_bins, energy_bins, photon_flux_density_bins);

	cout << "> Integration of test function yields " << result << endl;
	double exact_result = 0.;

	for(unsigned int i = 0; i < gamma0_list.size(); ++i){
		exact_result += beamParams[i]*(
				sqrt(PI/2.)*gamma0_list[i]*erf((energy_bins[NBINS - 1] - e0_list[2*i])/(sqrt(2)*gamma0_list[i])) - 
				sqrt(PI/2.)*gamma0_list[i]*erf((energy_bins[0] - e0_list[2*i])/(sqrt(2)*gamma0_list[i]))) * 
				(
			 	sqrt(PI/2.)*gamma0_list[i]*erf((energy_bins[NBINS - 1] - e0_list[2*i + 1])/(sqrt(2)*gamma0_list[i])) - 
				sqrt(PI/2.)*gamma0_list[i]*erf((energy_bins[0] - e0_list[2*i + 1])/(sqrt(2)*gamma0_list[i])));
	}	

	cout << "\tExact result: " << exact_result << " (" << ((result - exact_result)/exact_result*100.) << " %)" << endl;
}

void Target::calculateAbsorption(double (&energy_bins)[NBINS]){

	double absorption = integrateEZHistogram(energy_bins, z_bins, resonance_absorption_density_bins);

	cout << "TARGET #" << target_number << ": '" << target_name << "'" << endl;
	cout << "INT ALPHA dE dZ = " << absorption << endl;
}

void Target::plotVelocityDistribution(){

	stringstream filename;
	filename << target_name << "_velocity_distribution.pdf";
	stringstream canvasname;
	canvasname << target_name << "_canvas";
	TCanvas *canvas = new TCanvas(canvasname.str().c_str(), target_name.c_str(), 0, 0, 800, 500);
	TLegend *legend = new TLegend(V_PLOT_LEGEND_X1, V_PLOT_LEGEND_Y1, V_PLOT_LEGEND_X2, V_PLOT_LEGEND_Y2);

	// It is enough to plot the velocity distribution in the binning of one of the excited states, because it is the same for all states of a target.
	crossSection->plot_vdist(velocity_bins[0], vdist_bins[0], target_name, canvas, legend, "Velocity Distribution");

	legend->Draw();
	canvas->SaveAs(filename.str().c_str());
	delete canvas;
}

void Target::plotDopplerShift(double (&energy_bins)[NBINS]){

	stringstream filename;
	filename << target_name << "_doppler_shift.pdf";
	stringstream canvasname;
	canvasname << target_name << "_canvas";
	TCanvas *canvas = new TCanvas(canvasname.str().c_str(), target_name.c_str(), 0, 0, 800, 500);
	TLegend *legend = new TLegend(CS_PLOT_LEGEND_X1, CS_PLOT_LEGEND_Y1, CS_PLOT_LEGEND_X2, CS_PLOT_LEGEND_Y2);

	crossSection->plot_dopplershift(energy_bins, crosssection_bins, dopplercs_bins, target_name, canvas, legend, "Doppler-shifted cross section");

	legend->Draw();
	canvas->SaveAs(filename.str().c_str());
	delete canvas;
	
}

void Target::plotMassAttenuation(double (&energy_bins)[NBINS]){

	stringstream filename;
	filename << target_name << "_mass_attenuation.pdf";
	stringstream canvasname;
	canvasname << target_name << "_canvas";
	TCanvas *canvas = new TCanvas(canvasname.str().c_str(), target_name.c_str(), 0, 0, 800, 500);
	TLegend *legend = new TLegend(MU_PLOT_LEGEND_X1, MU_PLOT_LEGEND_Y1, MU_PLOT_LEGEND_X2, MU_PLOT_LEGEND_Y2);

	absorption->plot_massattenuation(energy_bins, massattenuation_bins, target_name, canvas, legend, "Mass attenuation coefficient");

	legend->Draw();
	canvas->SaveAs(filename.str().c_str());
	delete canvas;
	
}

void Target::plotMu(){

	stringstream filename;
	filename << target_name << "_mu.pdf";
	stringstream canvasname;
	canvasname << target_name << "_canvas";
	TCanvas *canvas = new TCanvas(canvasname.str().c_str(), target_name.c_str(), 0, 0, 800, 500);
	TLegend *legend = new TLegend(MU_PLOT_LEGEND_X1, MU_PLOT_LEGEND_Y1, MU_PLOT_LEGEND_X2, MU_PLOT_LEGEND_Y2);

	canvas->SetLogx();
	canvas->SetLogy();

	absorption->plot_total_massattenuation(target_name, canvas, legend, "Mass attenuation #mu");

	legend->Draw();
	canvas->SaveAs(filename.str().c_str());
	delete canvas;
	
}

void Target::plotPhotonFluxDensity(double (&energy_bins)[NBINS]){

	stringstream filename;
	filename << target_name << "_phi.pdf";
	stringstream canvasname;
	canvasname << target_name << "_canvas";
	TCanvas *canvas = new TCanvas(canvasname.str().c_str(), target_name.c_str(), 0, 0, 800, 500);
	TLegend *legend = new TLegend(MU_PLOT_LEGEND_X1, MU_PLOT_LEGEND_Y1, MU_PLOT_LEGEND_X2, MU_PLOT_LEGEND_Y2);

	absorption->plot_photon_flux_density(energy_bins, z_bins, photon_flux_density_bins, target_name, canvas, legend, "Photon flux density #Phi");

	legend->Draw();
	canvas->SaveAs(filename.str().c_str());
	delete canvas;
}

void Target::plotTestIntegration(double (&energy_bins)[NBINS]){

	stringstream filename;
	filename << target_name << "_function.pdf";
	stringstream canvasname;
	canvasname << target_name << "_canvas";
	TCanvas *canvas = new TCanvas(canvasname.str().c_str(), target_name.c_str(), 0, 0, 800, 500);
	TLegend *legend = new TLegend(MU_PLOT_LEGEND_X1, MU_PLOT_LEGEND_Y1, MU_PLOT_LEGEND_X2, MU_PLOT_LEGEND_Y2);

	absorption->plot_photon_flux_density(energy_bins, energy_bins, photon_flux_density_bins, target_name, canvas, legend, "Function inside the integral, ignore axis labels");

	legend->Draw();
	canvas->SaveAs(filename.str().c_str());
	delete canvas;
}

void Target::plotResonanceAbsorptionDensity(double (&energy_bins)[NBINS]){

	stringstream filename;
	filename << target_name << "_alpha.pdf";
	stringstream canvasname;
	canvasname << target_name << "_canvas";
	TCanvas *canvas = new TCanvas(canvasname.str().c_str(), target_name.c_str(), 0, 0, 800, 500);
	TLegend *legend = new TLegend(MU_PLOT_LEGEND_X1, MU_PLOT_LEGEND_Y1, MU_PLOT_LEGEND_X2, MU_PLOT_LEGEND_Y2);

	absorption->plot_resonance_absorption_density(energy_bins, z_bins, resonance_absorption_density_bins, target_name, canvas, legend, "Resonance absorption density #alpha");

	legend->Draw();
	canvas->SaveAs(filename.str().c_str());
	delete canvas;
	
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


double Target::normalizeVDist(unsigned int i){

	double sum = 0.;

	for(unsigned int j = 0; j < NBINS-1; j++){
		sum += vdist_bins[i][j]*(velocity_bins[i][j + 1] - velocity_bins[i][j]);
	}

	return 1./sum;
}
