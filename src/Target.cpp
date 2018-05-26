/*    
    This file is part of SeAN.

    SeAN is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SeAN is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SeAN.  If not, see <http://www.gnu.org/licenses/>.
*/


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

void Target::initialize(const vector<double> &energy_bins){
	unsigned int nenergies = (unsigned int) settings.energy[target_number].size();

	// Reserve space for the histograms
	for(unsigned int i = 0; i < nenergies; ++i){
		crosssection_at_rest_histogram.push_back(vector<double>(settings.nbins_e, 0.));
		velocity_distribution_bins.push_back(vector<double>(settings.nbins_e, 0.));
		velocity_distribution_histogram.push_back(vector<double>(settings.nbins_e, 0.));
	}

	crosssection_histogram = vector<double>(settings.nbins_e, 0.);
	incident_beam_histogram = vector<double>(settings.nbins_e, 0.);
	mass_attenuation_histogram = vector<double>(settings.nbins_e, 0.);
	mass_attenuation_file.push_back(vector<double>());
	mass_attenuation_file.push_back(vector<double>());

	z_bins = vector<double>(settings.nbins_z, 0.);

	for(unsigned int i = 0; i < settings.nbins_z; ++i){
		photon_flux_density_histogram.push_back(vector<double>(settings.nbins_e, 0.));
		resonance_absorption_density_histogram.push_back(vector<double>(settings.nbins_e, 0.));
	}

	energy_boosted.reserve(nenergies);

	// Shift resonance energies due to target velocity
	boost_and_recoil();
	// Calculate z bins
	calculateZBins();
	// Calculate cross section
	calculateCrossSectionAtRest(energy_bins);
	// Calculate velocity distribution
	calculateVelocityDistribution(energy_bins);
	// Calculate mass attenuation
	calculateMassAttenuation(energy_bins);

}

void Target::calculateCrossSection(const vector<double> &energy_bins){

	// Need to treat the four special cases
	// 1. Velocity distribution at T = 0 K
	// 2. Maxwell-Boltzmann with approximation
	// 3. Arbitrary doppler-shifted cross section or velocity distribution
	// 4. Phonon density distribution
	// separately, because it is possible to calculate the doppler-shifted cross section without the expensive double integral 
	if(settings.dopplerBroadening[target_number] == dopplerModel::zero){
		crossSection.no_dopplershift(crosssection_at_rest_histogram, crosssection_histogram);
	}

	else if(settings.dopplerBroadening[target_number] == dopplerModel::mba){
		crossSection.maxwell_boltzmann_approximation(energy_bins, crosssection_histogram, energy_boosted, target_number);
	}

	else if(settings.dopplerBroadening[target_number] == dopplerModel::mbad){
		crossSection.maxwell_boltzmann_approximation_debye(energy_bins, crosssection_histogram, energy_boosted, target_number);
	}

	else if(settings.dopplerBroadening[target_number] == dopplerModel::arb_cs){
		crossSection.arbitrary_cross_section(energy_bins, crosssection_histogram, energy_bins_file, cross_section_file);
	}
	else if(settings.dopplerBroadening[target_number] == dopplerModel::phdos){
		phononDensity.calculateCrossSection(energy_bins, energy_boosted, crosssection_histogram, omega_s_file, target_number);
	}

	else if(settings.exact){
		crossSection.integration_input(crosssection_at_rest_histogram, velocity_distribution_histogram);
		crossSection.dopplershift(energy_bins, crosssection_histogram, crosssection_at_rest_histogram, velocity_distribution_bins, velocity_distribution_histogram, vdist_norm, energy_boosted);
	} else{
		crossSection.fft_input(energy_bins, crosssection_at_rest_histogram, velocity_distribution_histogram, energy_boosted);
		crossSection.dopplershiftFFT(energy_bins, crosssection_histogram, crosssection_at_rest_histogram, velocity_distribution_bins, velocity_distribution_histogram, vdist_norm, vdist_centroid);
	}

	if(settings.verbosity > 0){
		crossSection.check_crosssection_normalization(energy_bins, crosssection_histogram, energy_boosted, target_number);
	}
}

void Target::boost_and_recoil(){
	double recoil = 0.;

	double beta = settings.velocity[target_number]/SPEEDOFLIGHT;
	
	unsigned int nenergies = (unsigned int) settings.energy[target_number].size();

	for(unsigned int i = 0; i < nenergies; ++i){
		if(settings.recoil){
			recoil = settings.energy[target_number][i]*settings.energy[target_number][i]/(2.*settings.mass[target_number]*AtomicMassUnit);

			cout << "RECOIL_CORRECTION:\t" << recoil << " eV" << endl;
		}

		energy_boosted.push_back((1. + beta)/sqrt(1. - beta*beta)*settings.energy[target_number][i] + recoil);
	}
}

void Target::calculateCrossSectionAtRest(const vector<double> &energy_bins){
	crossSection.breit_wigner(energy_bins, crosssection_at_rest_histogram, energy_boosted, target_number);
}

void Target::calculateVelocityDistribution(const vector<double> &energy_bins){

	if(settings.dopplerBroadening[target_number] != dopplerModel::arb_cs){
		crossSection.calculateVelocityBins(energy_bins, velocity_distribution_bins, energy_boosted, target_number);
	}

	stringstream filename;

	switch(settings.dopplerBroadening[target_number]){
		case dopplerModel::zero:
			crossSection.absolute_zero(velocity_distribution_bins, velocity_distribution_histogram, target_number);
			break;

		case dopplerModel::mb:
		case dopplerModel::mba:
		case dopplerModel::mbad:
			crossSection.maxwell_boltzmann(velocity_distribution_bins, velocity_distribution_histogram, target_number);
			break;

		case dopplerModel::mbd:
			crossSection.maxwell_boltzmann_debye(velocity_distribution_bins, velocity_distribution_histogram, target_number);
			break;

		case dopplerModel::arb_vdist:
			filename << VELOCITY_DISTRIBUTION_DIR << settings.velocityBinFile[target_number];
			inputReader.read1ColumnFile(velocity_bins_file, filename.str());

			filename.str("");
			filename.clear();
			filename << VELOCITY_DISTRIBUTION_DIR << settings.vDistFile[target_number];
			inputReader.read1ColumnFile(velocity_distribution_file, filename.str());

			crossSection.arbitrary_velocity_distribution(velocity_distribution_bins, velocity_distribution_histogram, velocity_bins_file, velocity_distribution_file, energy_boosted, target_number);
			break;

		case dopplerModel::arb_cs:
			filename << CROSS_SECTION_DIR << settings.energyBinFile[target_number];
			inputReader.read1ColumnFile(energy_bins_file, filename.str());

			filename.str("");
			filename.clear();
			filename << CROSS_SECTION_DIR << settings.crosssectionFile[target_number];
			inputReader.read1ColumnFile(cross_section_file, filename.str());
			break;

		case dopplerModel::phdos:
			filename << PHONON_DIR << settings.omegaFile[target_number];
			inputReader.read1ColumnFile(omega_s_file, filename.str());

		//	filename.str("");
		//	filename.clear();
		//	filename << PHONON_DIR << settings.polarizationFile[target_number];
		//	inputReader.read3ColumnFile(e_s_file, filename.str());

		//	filename.str("");
		//	filename.clear();
		//	filename << PHONON_DIR << settings.momentumFile[target_number];
		//	inputReader.read3ColumnFile(p_file, filename.str());
	}

	vDistInfo();
}

void Target::plot(vector<double> &energy_bins, unsigned int n_setting) {

	stringstream filename;
	
	filename << settings.targetNames[target_number] << "_crosssection_at_rest_" << n_setting;

	plotter.plotMultiple1DHistogramsAndSum(energy_bins, crosssection_at_rest_histogram, filename.str(), "Energy / eV", "Cross section / fm^{2}");

	// Plot velocity distribution
	// In fact, each resonance has its own binning, but plot only the velocity distribution for the first one since they only differ in the binning
	filename.str("");
	filename.clear();

	filename << settings.targetNames[target_number] << "_velocity_distribution_" << n_setting;

	plotter.plot1DHistogram(velocity_distribution_bins[0], velocity_distribution_histogram[0], filename.str(), "Velocity / c", "Velocity distribution");

	// Plot cross section
	filename.str("");
	filename.clear();

	filename << settings.targetNames[target_number] << "_crosssection_" << n_setting;
	plotter.plot1DHistogram(energy_bins, crosssection_histogram, filename.str(), "Energy / eV", "Cross section / fm^{2}");

	// Plot incident beam
	filename.str("");
	filename.clear();

	filename << settings.targetNames[target_number] << "_incident_beam_" << n_setting;
	plotter.plot1DHistogram(energy_bins, incident_beam_histogram, filename.str(), "Energy / eV", "Beam intensity / a.u.");

	// Plot mass attenuation
	filename.str("");
	filename.clear();

	filename << settings.targetNames[target_number] << "_mass_attenuation_" << n_setting;
	plotter.plot1DHistogram(energy_bins, mass_attenuation_histogram, filename.str(), "Energy / eV", "Mass attenuation / fm^{2} / atom");

	// Plot photon flux density
	filename.str("");
	filename.clear();

	filename << settings.targetNames[target_number] << "_photon_flux_density_" << n_setting;
	plotter.plot2DHistogram(z_bins, energy_bins, photon_flux_density_histogram, filename.str(), "z / atoms/fm^{2}", "Energy / eV", "#phi");
	
	// Plot resonance absorptions density
	filename.str("");
	filename.clear();

	filename << settings.targetNames[target_number] << "_resonance_absorption_density_" << n_setting;
	plotter.plot2DHistogram(z_bins, energy_bins, resonance_absorption_density_histogram, filename.str(), "z / atoms/fm^{2}", "Energy / eV", "#alpha / fm^{2}");
}

void Target::calculateIncidentBeam(const vector<double> &energy_bins){
	
	switch(settings.incidentBeam){
		case incidentBeamModel::constant:
			absorption.const_beam(energy_bins, incident_beam_histogram);
			break;
		case incidentBeamModel::gauss:
			absorption.gauss_beam(energy_bins, incident_beam_histogram);
			break;
		case incidentBeamModel::arb:
			stringstream filename;
			filename << BEAM_DIR << settings.incidentBeamFile;
			inputReader.read2ColumnFile(incident_beam_file, filename.str());
			absorption.arbitrary_beam(energy_bins, incident_beam_histogram, incident_beam_file);
			break;
	}
};

void Target::calculateIncidentBeam(const vector< vector<double> > &photon_flux_density_histogram){

	for(unsigned int i = 0; i < settings.nbins_e; ++i){
		incident_beam_histogram[i] = photon_flux_density_histogram[settings.nbins_z - 1][i];
	}
}

void Target::calculateZBins(){
	
	double delta_z = settings.thickness[target_number]/settings.nbins_z;

	for(unsigned int i = 0; i < settings.nbins_z; ++i){
		z_bins[i] = i*delta_z;
	}
}

void Target::calculateMassAttenuation(const vector<double> &energy_bins){
	
	stringstream filename;

	switch(settings.mAtt[target_number]){
		case mAttModel::constant:
			absorption.const_mass_attenuation(mass_attenuation_histogram, target_number);
			break;
		case mAttModel::arb:
			filename << MU_DIR << settings.mAttFile[target_number];
			inputReader.read2ColumnFile(mass_attenuation_file, filename.str()); 
			absorption.arbitrary_mass_attenuation(energy_bins, mass_attenuation_file, mass_attenuation_histogram);
			break;
		case mAttModel::nist:
			filename << settings.mAttFile[target_number];
			inputReader.readNIST(mass_attenuation_file, filename.str()); 
			absorption.nist_mass_attenuation(energy_bins, mass_attenuation_file, mass_attenuation_histogram, target_number);
			break;
	}
}

void Target::calculateTransmission(const vector<double> energy_bins){
	
	absorption.photon_flux_density(crosssection_histogram, mass_attenuation_histogram, z_bins, incident_beam_histogram, photon_flux_density_histogram);

	absorption.resonance_absorption_density(crosssection_histogram, photon_flux_density_histogram, resonance_absorption_density_histogram);
}

void Target::write(const vector<double> &energy_bins, const unsigned int n_setting) {
	
	stringstream filename;

	// Write energy bins
	filename << settings.targetNames[target_number] << "_energy_bins_" << n_setting;
	writer.write1DHistogram(energy_bins, filename.str(), "Energy / eV");
	
	// Write cross section at rest
	
	for(unsigned int i = 0; i < energy_boosted.size(); ++i){
		filename.str("");
		filename.clear();

		filename << settings.targetNames[target_number] << "_crosssection_at_rest_"  << n_setting << "_" << i;

		writer.write1DHistogram(crosssection_at_rest_histogram[i], filename.str(), "Cross section / fm^2");
		writer.write1DCalibration(energy_bins, CAL_FILE_NAME, filename.str());
	}

	// Write cross section
	filename.str("");
	filename.clear();

	filename << settings.targetNames[target_number] << "_crosssection_" << n_setting;

	writer.write1DHistogram(crosssection_histogram, filename.str(), "Cross section / fm^2");
	writer.write1DCalibration(energy_bins, CAL_FILE_NAME, filename.str());

	// Write velocity distribution 
	// In fact, each resonance has its own binning, but write only the velocity distribution for the first one
	
	filename.str("");
	filename.clear();

	filename << settings.targetNames[target_number] << "_velocity_bins_" << n_setting;
	writer.write1DHistogram(velocity_distribution_bins[0], filename.str(), "Velocity / c");
	writer.write1DCalibration(velocity_distribution_bins[0], CAL_FILE_NAME, filename.str());

	filename.str("");
	filename.clear();

	filename << settings.targetNames[target_number] << "_velocity_histogram_" << n_setting;
	writer.write1DHistogram(velocity_distribution_histogram[0], filename.str(), "Velocity distribution");

	// Write incident beam

	filename.str("");
	filename.clear();

	filename << settings.targetNames[target_number] << "_incident_beam_" << n_setting;
	writer.write1DHistogram(incident_beam_histogram, filename.str(), "Beam intensity distribution / a.u.");
	writer.write1DCalibration(energy_bins, CAL_FILE_NAME, filename.str());

	// Write mass attenuation
	filename.str("");
	filename.clear();

	filename << settings.targetNames[target_number] << "_mass_attenuation_" << n_setting;
	writer.write1DHistogram(mass_attenuation_histogram, filename.str(), "Mass attenuation / fm^2 / atom");
	writer.write1DCalibration(energy_bins, CAL_FILE_NAME, filename.str());
	
	// Write photon flux density
	filename.str("");
	filename.clear();

	filename << settings.targetNames[target_number] << "_photon_flux_density_" << n_setting;
	writer.write2DHistogram(photon_flux_density_histogram, filename.str(), "Phi (z, E = const)", "Phi (z = const, E)");
	
	// Write resonance absorption density
	filename.str("");
	filename.clear();

	filename << settings.targetNames[target_number] << "_resonance_absorption_density_" << n_setting;
	writer.write2DHistogram(resonance_absorption_density_histogram, filename.str(), "Alpha (z, E = const) / fm^2", "Alpha (z = const, E) / fm^2");
}

void Target::write_results(string outputfile) const{

	stringstream filename;
	filename << "output/" << outputfile;

	ofstream ofile;
	ofile.open(filename.str(), std::ios_base::out | std::ios_base::app);

        if(!ofile.is_open()){
                cout << "Error: " << __FILE__ << ":" << __LINE__ << ": "; 
		cout << " write_results(): File '" << outputfile << "' could not be opened." << endl;
		abort();
	}

	ofile << "Target #" << target_number << "\t" << n_resonantly_scattered << endl;
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

void Target::calculateResonantScattering(const vector<double> energy_bins){

	n_resonantly_scattered = integrator.integrate2DHistogram(z_bins, energy_bins, resonance_absorption_density_histogram);

}

void Target::print_results(){
	cout << settings.targetNames[target_number] << "\t" << n_resonantly_scattered << endl;
}
