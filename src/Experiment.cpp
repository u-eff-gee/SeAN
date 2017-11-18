#include <iostream>
#include <fstream>
#include <string>
#include <regex>

#include "Experiment.h"

using std::cout;
using std::endl;
using std::ifstream;
using std::string;
using std::regex;
using std::regex_replace;

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
		targets.push_back(new Target(i, settings));
		targets[i]->initialize(energy_bins);
	}
}

void Experiment::crossSections(){

	unsigned int ntargets = (unsigned int) settings.targetNames.size();	

	for(unsigned int i = 0; i < ntargets; ++i){
		targets[i]->calculateCrossSection(energy_bins);
	}
};

void Experiment::transmission(){

	unsigned int ntargets = (unsigned int) settings.targetNames.size();	

	for(unsigned int i = 0; i < ntargets; ++i){
		// The incident beam on the first target is user-defined, the rest is determined by the transmission of the previous target
		if(i == 0){
			targets[i]->calculateIncidentBeam(energy_bins);
		} else{
			;
		}
	}
	//targets[0]->calculateIncidentBeam(energy_bins, beam_ID, beamParams);
//	for(unsigned int i = 0; i < targets.size(); ++i){
//		targets[i]->calculateMassAttenuation(energy_bins);
//		targets[i]->calculateZBins();
//		if(i > 0){
//			targets[i]->setIncidentBeam(targets[i-1]->getTransmittedBeam());
//		}
//		targets[i]->calculatePhotonFluxDensity();
//		targets[i]->calculateTransmittedBeam();
//		targets[i]->calculateResonanceAbsorptionDensity();
//		targets[i]->calculateAbsorption(energy_bins);
//
//		if(settings.plot){
//			targets[i]->plotMassAttenuation(energy_bins);
//			targets[i]->plotMu();
//			targets[i]->plotPhotonFluxDensity(energy_bins);
//			targets[i]->plotResonanceAbsorptionDensity(energy_bins);
//		}
//		if(settings.write){
//			targets[i]->write(energy_bins);
//		}
//	}
}

void Experiment::plot(){

	unsigned int ntargets = (unsigned int) settings.targetNames.size();	

	for(unsigned int i = 0; i < ntargets; ++i){
		targets[i]->plot(energy_bins);
	}
}

void Experiment::write(){

	unsigned int ntargets = (unsigned int) settings.targetNames.size();	

	for(unsigned int i = 0; i < ntargets; ++i){
		targets[i]->write(energy_bins);
	}
}
