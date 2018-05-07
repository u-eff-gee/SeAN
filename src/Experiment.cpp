#include <iostream>
#include <fstream>
#include <string>
#include <regex>
#include <sstream>

#include "Experiment.h"
#include "Config.h"

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::string;
using std::regex;
using std::regex_replace;
using std::stringstream;

Experiment::Experiment(Settings &s){
	settings = s;
}

void Experiment::initialize(){
	createEnergyBins(settings.emin, settings.emax);
	createTargets();
}

void Experiment::createEnergyBins(double emin, double emax){

	energy_bins.reserve(settings.nbins_e);

	double delta_e = (double) (emax - emin)/(settings.nbins_e - 1);

	for(unsigned int i = 0; i < settings.nbins_e; ++i)
		energy_bins.push_back(emin + i*delta_e);
};

void Experiment::createTargets(){

	unsigned int ntargets = (unsigned int) settings.targetNames.size();	

	for(unsigned int i = 0; i < ntargets; ++i){
		targets.push_back(Target(i, settings));
		targets[i].initialize(energy_bins);
	}
}

void Experiment::crossSections(){

	unsigned int ntargets = (unsigned int) settings.targetNames.size();	

	for(unsigned int i = 0; i < ntargets; ++i){
		targets[i].calculateCrossSection(energy_bins);
	}
};

void Experiment::transmission(){

	unsigned int ntargets = (unsigned int) settings.targetNames.size();	

	for(unsigned int i = 0; i < ntargets; ++i){
		// The incident beam on the first target is user-defined, the rest is determined by the transmission of the previous target
		if(i == 0){
			targets[0].calculateIncidentBeam(energy_bins);
			targets[0].calculateTransmission(energy_bins);
		} else{
			targets[i].calculateIncidentBeam(targets[i-1].getPhotonFluxDensity());
			targets[i].calculateTransmission(energy_bins);
		}
	}
}

void Experiment::resonant_scattering(){

	unsigned int ntargets = (unsigned int) settings.targetNames.size();	

	for(unsigned int i = 0; i < ntargets; ++i){
		targets[i].calculateResonantScattering(energy_bins);
	}
}

void Experiment::print_results(){

	unsigned int ntargets = (unsigned int) settings.targetNames.size();	

	cout << HORIZONTAL_LINE << endl;
	cout << ">>> SeAN RESULTS" << endl;
	cout << "TARGET NAME\tRESONANT SCATTERING" << endl;

	for(unsigned int i = 0; i < ntargets; ++i){
		targets[i].print_results();
	}
}

void Experiment::plot(unsigned int n_setting) {

	unsigned int ntargets = (unsigned int) settings.targetNames.size();	

	for(unsigned int i = 0; i < ntargets; ++i){
		targets[i].plot(energy_bins, n_setting);
	}
}

void Experiment::write(unsigned int n_setting) {

	unsigned int ntargets = (unsigned int) settings.targetNames.size();	

	for(unsigned int i = 0; i < ntargets; ++i){
		targets[i].write(energy_bins, n_setting);
	}
}

void Experiment::write_results(string outputfile, unsigned int n_setting) const{

	stringstream filename;
	filename << "output/" << outputfile;

	ofstream ofile;
	ofile.open(filename.str(), std::ios_base::out | std::ios_base::app);

        if(!ofile.is_open()){
                cout << "Error: " << __FILE__ << ":" << __LINE__ << ": "; 
		cout << " write_results(): File '" << outputfile << "' could not be opened." << endl;
		abort();
	}
	
	unsigned int ntargets = (unsigned int) settings.targetNames.size();	

	ofile << ">>> SeAN RESULTS #" << n_setting << endl;
	ofile << HORIZONTAL_LINE << endl;
	ofile << "TARGET NAME\tRESONANT SCATTERING" << endl;

	for(unsigned int i = 0; i < ntargets; ++i){
		targets[i].write_results(outputfile);
	}
}
