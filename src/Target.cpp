#include "Target.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>

#include "TCanvas.h"
#include "TLegend.h"

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::stringstream;
using std::regex;
using std::scientific;

Target::~Target(){
	delete crossSection;
	delete absorption;

	delete &crosssection_at_rest_histogram;
	delete &velocity_distribution_histogram;
	delete &velocity_distribution_histogram;
	delete &z_bins;
	delete &energy_boosted;
	delete &vdist_norm;
	delete &vdist_centroid;

	delete &photon_flux_density_histogram;
	delete &resonance_absorption_density_histogram;

	delete &incident_beam_histogram;
	delete &crosssection_histogram;
	delete &massattenuation_histogram;
	delete &transmitted_beam_histogram;
};

void Target::initialize(vector<double> &energy_bins){
	unsigned int nenergies = (unsigned int) settings.energy[target_number].size();

	// Reserve space for the histograms
	for(unsigned int i = 0; i < nenergies; ++i){
		crosssection_at_rest_histogram.push_back(vector<double>(settings.nbins_e, 0.));
		velocity_distribution_bins.push_back(vector<double>(settings.nbins_e, 0.));
		velocity_distribution_histogram.push_back(vector<double>(settings.nbins_e, 0.));
	}

	crosssection_histogram = vector<double>(settings.nbins_e, 0.);
	incident_beam_histogram = vector<double>(settings.nbins_e, 0.);

	z_bins = vector<double>(settings.nbins_z, 0.);

	photon_flux_density_histogram.reserve(settings.nbins_e*settings.nbins_z);
	resonance_absorption_density_histogram.reserve(settings.nbins_e*settings.nbins_z);

	energy_boosted.reserve(nenergies);

	// Initialize input reader and (if necessary) read mass attenuation
	inputReader = new InputReader();

	// Initialize calculators
	crossSection = new CrossSection(settings);
	absorption = new Absorption(settings);

	// Initialize plotter
	plotter = new Plotter();

	// Initialize writer
	writer = new Writer();

	// Shift resonance energies due to target velocity
	boostEnergies();
	// Calculate z bins
	calculateZBins();
	// Calculate cross section
	calculateCrossSectionAtRest(energy_bins);
	// Calculate velocity distribution
	calculateVelocityDistribution(energy_bins);

}

void Target::calculateCrossSection(vector<double> &energy_bins){

	// Need to treat the two special cases
	// 1. Velocity distribution at T = 0 K
	// 2. Maxwell-Boltzmann with approximation
	// separately, because it is possible to calculate the doppler-shifted cross section without the expensive double integral 
	if(settings.vDist[target_number] == vDistModel::zero){
		;
	}

	else if(settings.vDist[target_number] == vDistModel::mba){
		;
	}

	else if(settings.exact){
		crossSection->integration_input(crosssection_at_rest_histogram, velocity_distribution_histogram);
		crossSection->dopplershift(energy_bins, crosssection_histogram, crosssection_at_rest_histogram, velocity_distribution_bins, velocity_distribution_histogram, vdist_norm, energy_boosted);
	} else{
		crossSection->fft_input(energy_bins, crosssection_at_rest_histogram, velocity_distribution_histogram, energy_boosted);
		crossSection->dopplershiftFFT(energy_bins, crosssection_histogram, crosssection_at_rest_histogram, velocity_distribution_bins, velocity_distribution_histogram, vdist_norm, vdist_centroid);
	}
}

void Target::boostEnergies(){
	double beta = settings.velocity[target_number]/SPEEDOFLIGHT;
	
	unsigned int nenergies = (unsigned int) settings.energy[target_number].size();

	for(unsigned int i = 0; i < nenergies; ++i)
		energy_boosted.push_back((1. + beta)/sqrt(1. - beta*beta)*settings.energy[target_number][i]);
}

void Target::calculateCrossSectionAtRest(vector<double> &energy_bins){
	crossSection->breit_wigner(energy_bins, crosssection_at_rest_histogram, energy_boosted, target_number);
}

void Target::calculateVelocityDistribution(vector<double> &energy_bins){

	crossSection->calculateVelocityBins(energy_bins, velocity_distribution_bins, energy_boosted, target_number);

	switch(settings.vDist[target_number]){
		case vDistModel::zero:
			for(unsigned int i = 0; i < crosssection_at_rest_histogram.size(); ++i){
				for(unsigned int j = 0; j < settings.nbins_e; ++j){
					crosssection_histogram[i] += crosssection_at_rest_histogram[i][j];
				}
			}
			break;

		case vDistModel::mb:
		case vDistModel::mba:
			crossSection->maxwell_boltzmann(velocity_distribution_bins, velocity_distribution_histogram, target_number);
			break;

		case vDistModel::arb:
			stringstream filename;
			filename << "velocity_distribution/" << settings.vDistFile[target_number];
			inputReader->read2ColumnFile(velocity_distribution_file, filename.str());
			crossSection->arbitrary_velocity_distribution(velocity_distribution_bins, velocity_distribution_histogram, velocity_distribution_file, energy_boosted, target_number);
			break;
	}

	vDistInfo();
}

void Target::plot(vector<double> &energy_bins){

	stringstream filename;
	
	filename << settings.targetNames[target_number] << "_crosssection_at_rest";

	plotter->plotMultiple1DHistogramsAndSum(energy_bins, crosssection_at_rest_histogram, filename.str());

	// Plot velocity distribution
	// In fact, each resonance has its own binning, but plot only the velocity distribution for the first one since they only differ in the binning
	filename.str("");
	filename.clear();

	filename << settings.targetNames[target_number] << "_velocity_distribution";

	plotter->plot1DHistogram(velocity_distribution_bins[0], velocity_distribution_histogram[0], filename.str());

	// Plot cross section
	filename.str("");
	filename.clear();

	filename << settings.targetNames[target_number] << "_crosssection";
	plotter->plot1DHistogram(energy_bins, crosssection_histogram, filename.str());

	// Plot incident beam
	filename.str("");
	filename.clear();

	filename << settings.targetNames[target_number] << "_incident_beam";
	plotter->plot1DHistogram(energy_bins, incident_beam_histogram, filename.str());
}

void Target::calculateIncidentBeam(const vector<double> &energy_bins){
	
	switch(settings.incidentBeam){
		case incidentBeamModel::constant:
			absorption->const_beam(energy_bins, incident_beam_histogram);
			break;
		case incidentBeamModel::gauss:
			absorption->gauss_beam(energy_bins, incident_beam_histogram);
			break;
		case incidentBeamModel::arb:
			stringstream filename;
			filename << "beam/" << settings.incidentBeamFile;
			inputReader->read2ColumnFile(incident_beam_file, filename.str());
			absorption->arbitrary_beam(energy_bins, incident_beam_histogram, incident_beam_file);
			break;
	}
};

void Target::calculateZBins(){
	
	double delta_z = settings.thickness[target_number]/settings.nbins_z;

	for(unsigned int i = 0; i < settings.nbins_z; ++i){
		z_bins[i] = i*delta_z;
	}
}

//void Target::calculateMassAttenuation(vector<double> &energy_bins){
//	if(massAttenuation_ID == "0"){
//	;} else{
//		//absorption->read_massattenuation_NIST(energy_bins, massattenuation_bins, massAttenuation_ID, mass);
//	}
//}
//
//void Target::setIncidentBeam(double &trans_beam_bins){
//	for(unsigned int i = 0; i < settings.nbins_e; ++i)
//		incident_beam_bins[i] = (&trans_beam_bins)[i];
//}
//
//void Target::calculatePhotonFluxDensity(){
//	absorption->photon_flux_density(dopplercs_bins, massattenuation_bins, z_bins, incident_beam_bins, photon_flux_density_bins);
//}
//
//void Target::calculateTransmittedBeam(){
//	for(unsigned int i = 0; i < settings.nbins_e; ++i)
//		transmitted_beam_bins[i] = photon_flux_density_bins[i][NBINS_Z - 1];
//}
//
//void Target::calculateResonanceAbsorptionDensity(){
//	absorption->resonance_absorption_density(dopplercs_bins, photon_flux_density_bins, resonance_absorption_density_bins);
//}
//
//double Target::integrateEZHistogram(vector<double> &energy_bins, vector<double> &z_bins, vector<vector<double> > &ezhist){
//
//	// Area of a bin in 2D plane
//	double bin_area = (energy_bins[1] - energy_bins[0])*(z_bins[1] - z_bins[0]);
//
//	// Crude implementation of the Riemann integral as a sum of bin contents times their dimension in energy- and z-direction. Since there are only settings.nbins_e-1 spaces between settings.nbins_e bins, leave out the last bin in each loop.
//	double integral = 0.;
//
//	#pragma omp parallel for reduction (+:integral)
//	for(unsigned int i = 0; i < settings.nbins_e - 1; ++i){
//		for(unsigned int j = 0; j < NBINS_Z - 1; ++j){
//			integral += ezhist[i][j]; 
//		}
//	}
//
//	// Implementation that starts at the 1st and not the 0th bin to check the validity of this approximation
////	double integral = 0.;
////	for(int i = 1; i < settings.nbins_e; ++i){
////		for(int j = 1; j < NBINS_Z; ++j){
////			integral += bin_area*ezhist[i][j]; 
////		}
////	}
//
//	return bin_area*integral;
//}
//
//double Target::integrateEEHistogram(vector<double> &energy_bins, vector<vector<double> > &eehist){
//
//	// Area of a bin in 2D plane
//	double bin_area = (energy_bins[1] - energy_bins[0])*(energy_bins[1] - energy_bins[0]);
//
//	// Crude implementation of the Riemann integral as a sum of bin contents times their dimension in energy- and z-direction. Since there are only settings.nbins_e-1 spaces between settings.nbins_e bins, leave out the last bin in each loop.
//	double integral = 0.;
//
//	for(unsigned int i = 0; i < settings.nbins_e - 1; ++i){
//		for(unsigned int j = 0; j < settings.nbins_e - 1; ++j){
//			integral += bin_area*eehist[i][j]; 
//		}
//	}
//
//	return integral;
//}
//
//void Target::calculateAbsorption(vector<double> &energy_bins){
//
//	double absorption = integrateEZHistogram(energy_bins, z_bins, resonance_absorption_density_bins);
//
//	cout << "TARGET #" << target_number << ": '" << target_name << "'" << endl;
//	cout << "INT ALPHA dE dZ = " << absorption << endl;
//}
//
//void Target::plotVelocityDistribution(){
//
//	stringstream filename;
//	filename << PLOT_OUTPUT_DIR << target_name << "_velocity_distribution.pdf";
//	stringstream canvasname;
//	canvasname << target_name << "_canvas";
//	TCanvas *canvas = new TCanvas(canvasname.str().c_str(), target_name.c_str(), 0, 0, 800, 500);
//	TLegend *legend = new TLegend(V_PLOT_LEGEND_X1, V_PLOT_LEGEND_Y1, V_PLOT_LEGEND_X2, V_PLOT_LEGEND_Y2);
//
//	// It is enough to plot the velocity distribution in the binning of one of the excited states, because it is the same for all states of a target.
//	crossSection->plot_vdist(velocity_bins[0], vdist_bins[0], target_name, canvas, legend, "Velocity Distribution");
//
//	legend->Draw();
//	canvas->SaveAs(filename.str().c_str());
//	delete canvas;
//}
//
//void Target::plotDopplerShift(vector<double> &energy_bins){
//
//	stringstream filename;
//	filename << PLOT_OUTPUT_DIR << target_name << "_doppler_shift.pdf";
//	stringstream canvasname;
//	canvasname << target_name << "_canvas";
//	TCanvas *canvas = new TCanvas(canvasname.str().c_str(), target_name.c_str(), 0, 0, 800, 500);
//	TLegend *legend = new TLegend(CS_PLOT_LEGEND_X1, CS_PLOT_LEGEND_Y1, CS_PLOT_LEGEND_X2, CS_PLOT_LEGEND_Y2);
//
//	crossSection->plot_dopplershift(energy_bins, crosssection_bins, dopplercs_bins, target_name, canvas, legend, "Doppler-shifted cross section");
//
//	legend->Draw();
//	canvas->SaveAs(filename.str().c_str());
//	delete canvas;
//	
//}
//
//void Target::plotMassAttenuation(vector<double> &energy_bins){
//
//	stringstream filename;
//	filename << PLOT_OUTPUT_DIR << target_name << "_mass_attenuation.pdf";
//	stringstream canvasname;
//	canvasname << target_name << "_canvas";
//	TCanvas *canvas = new TCanvas(canvasname.str().c_str(), target_name.c_str(), 0, 0, 800, 500);
//	TLegend *legend = new TLegend(MU_PLOT_LEGEND_X1, MU_PLOT_LEGEND_Y1, MU_PLOT_LEGEND_X2, MU_PLOT_LEGEND_Y2);
//
//	absorption->plot_massattenuation(energy_bins, massattenuation_bins, target_name, canvas, legend, "Mass attenuation coefficient");
//
//	legend->Draw();
//	canvas->SaveAs(filename.str().c_str());
//	delete canvas;
//	
//}
//
//void Target::plotMu(){
//
//	stringstream filename;
//	filename << PLOT_OUTPUT_DIR << target_name << "_mu.pdf";
//	stringstream canvasname;
//	canvasname << target_name << "_canvas";
//	TCanvas *canvas = new TCanvas(canvasname.str().c_str(), target_name.c_str(), 0, 0, 800, 500);
//	TLegend *legend = new TLegend(MU_PLOT_LEGEND_X1, MU_PLOT_LEGEND_Y1, MU_PLOT_LEGEND_X2, MU_PLOT_LEGEND_Y2);
//
//	canvas->SetLogx();
//	canvas->SetLogy();
//
//	absorption->plot_total_massattenuation(target_name, canvas, legend, "Mass attenuation #mu");
//
//	legend->Draw();
//	canvas->SaveAs(filename.str().c_str());
//	delete canvas;
//	
//}
//
//void Target::plotPhotonFluxDensity(vector<double> &energy_bins){
//
//	stringstream filename;
//	filename << PLOT_OUTPUT_DIR << target_name << "_phi.pdf";
//	stringstream canvasname;
//	canvasname << target_name << "_canvas";
//	TCanvas *canvas = new TCanvas(canvasname.str().c_str(), target_name.c_str(), 0, 0, 800, 500);
//	TLegend *legend = new TLegend(MU_PLOT_LEGEND_X1, MU_PLOT_LEGEND_Y1, MU_PLOT_LEGEND_X2, MU_PLOT_LEGEND_Y2);
//
//	absorption->plot_photon_flux_density(energy_bins, z_bins, photon_flux_density_bins, target_name, canvas, legend, "Photon flux density #Phi");
//
//	legend->Draw();
//	canvas->SaveAs(filename.str().c_str());
//	delete canvas;
//}
//
//void Target::plotResonanceAbsorptionDensity(vector<double> &energy_bins){
//
//	stringstream filename;
//	filename << PLOT_OUTPUT_DIR << target_name << "_alpha.pdf";
//	stringstream canvasname;
//	canvasname << target_name << "_canvas";
//	TCanvas *canvas = new TCanvas(canvasname.str().c_str(), target_name.c_str(), 0, 0, 800, 500);
//	TLegend *legend = new TLegend(MU_PLOT_LEGEND_X1, MU_PLOT_LEGEND_Y1, MU_PLOT_LEGEND_X2, MU_PLOT_LEGEND_Y2);
//
//	absorption->plot_resonance_absorption_density(energy_bins, z_bins, resonance_absorption_density_bins, target_name, canvas, legend, "Resonance absorption density #alpha");
//
//	legend->Draw();
//	canvas->SaveAs(filename.str().c_str());
//	delete canvas;
//	
//}
//
void Target::write(vector<double> &energy_bins){
	
	stringstream filename;

	// Write energy bins
	filename << settings.targetNames[target_number] << "_energy_bins";
	writer->write1DHistogram(energy_bins, filename.str(), "Energy / eV");
	
	// Write cross section at rest
	
	for(unsigned int i = 0; i < energy_boosted.size(); ++i){
		filename.str("");
		filename.clear();

		filename << settings.targetNames[target_number] << "_crosssection_at_rest_" << i;

		writer->write1DHistogram(crosssection_at_rest_histogram[i], filename.str(), "Cross section / fm^2");
	}

	// Write velocity distribution 
	// In fact, each resonance has its own binning, but write only the velocity distribution for the first one
	
	filename.str("");
	filename.clear();

	filename << settings.targetNames[target_number] << "_velocity_bins";
	writer->write1DHistogram(velocity_distribution_bins[0], filename.str(), "Velocity / c");

	filename.str("");
	filename.clear();

	filename << settings.targetNames[target_number] << "_velocity_histogram";
	writer->write1DHistogram(velocity_distribution_histogram[0], filename.str(), "Velocity distribution");

	// Write incident beam

	filename.str("");
	filename.clear();

	filename << settings.targetNames[target_number] << "_incident_beam";
	writer->write1DHistogram(incident_beam_histogram, filename.str(), "Beam intensity distribution");

}

void Target::vDistInfo(){
	
//	 Determine the integral of the velocity distribution and the centroid of the distribution.
//	 The integral is needed to normalize the pseudo-convolution with the cross section, since this should conserve the area under the resonance curve.
//	 The centroid is needed for the FFT convolution. The algorithm which is used re-orders the input array by the shift of the centroid of the velocity distribution from the array centroid.
	double norm = 0.;
	double sum = 0.;
	double weightedsum = 0.;
	double velocity_bin_low = 0.;
	double velocity_bin_high = 0.;
	double vdist = 0.;

	for(unsigned int i = 0; i < velocity_distribution_bins.size(); ++i){

		sum = 0.;
		norm = 0.;
		weightedsum = 0.;

		for(unsigned int j = 0; j < settings.nbins_e - 1; ++j){
			vdist = velocity_distribution_histogram[i][j];
			velocity_bin_low = velocity_distribution_bins[i][j];
			velocity_bin_high = velocity_distribution_bins[i][j + 1];
			norm += vdist*(velocity_bin_high - velocity_bin_low);
			weightedsum += vdist*j;
			sum += vdist;
		}

		vdist_norm.push_back(fabs(norm));
		vdist_centroid.push_back((unsigned int) (weightedsum/sum));
	}
}
