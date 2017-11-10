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
	settings.nbins_e = 0;
	settings.nbins_z = 0;
}

void Experiment::initialize(){
	createEnergyBins(settings.emin, settings.emax);
	createTargets();
}

void Experiment::createEnergyBins(double emin, double emax){

	energy_bins.reserve(settings.nbins_e);

	double delta_e = (double) (emax - emin)/(NBINS - 1);

	for(unsigned int i = 0; i < NBINS; ++i)
		energy_bins.push_back(emin + i*delta_e);
};

void Experiment::createTargets(){

	unsigned int ntargets = (unsigned int) settings.targetNames.size();	

	for(unsigned int i = 0; i < ntargets; ++i){
		targets.push_back(new Target(i, settings));
		targets[i]->initialize();
	}
}

void Experiment::crossSections(){
//	for(unsigned int i = 0; i < targets.size(); ++i){
//		targets[i]->calculateCrossSection(energy_bins);
//		if(settings.exact){
//			targets[i]->calculateVelocityDistribution(energy_bins);
//			targets[i]->calculateDopplerShift(energy_bins);
//		} else{
//			targets[i]->calculateVelocityDistribution(energy_bins);
//			targets[i]->calculateDopplerShiftFFT(energy_bins);
//		}
//
//		if(settings.plot){
//			targets[i]->plotCrossSection(energy_bins);
//			targets[i]->plotVelocityDistribution();
//			targets[i]->plotDopplerShift(energy_bins);
//		}
//	}
};

void Experiment::transmission(){
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
