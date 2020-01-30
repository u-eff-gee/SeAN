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
		if(settings.direct){
			targets[i].calculateCrossSectionDirectly(energy_bins);
		} else{
			targets[i].calculateCrossSection(energy_bins);
		}
	}
};

void Experiment::transmission(){

	unsigned int ntargets = (unsigned int) settings.targetNames.size();	

	for(unsigned int i = 0; i < ntargets; ++i){
		// The incident beam on the first target is user-defined, the rest is determined by the transmission of the previous target
		if(i == 0){
			targets[0].calculateIncidentBeam(energy_bins);
			targets[0].calculateTransmission();
		} else{
			targets[i].calculateIncidentBeam(targets[i-1].getPhotonFluxDensity());
			targets[i].calculateTransmission();
		}
	}
}

void Experiment::resonant_scattering(){

	unsigned int ntargets = (unsigned int) settings.targetNames.size();	

	for(unsigned int i = 0; i < ntargets; ++i){
		targets[i].calculateResonantScattering(energy_bins);
		if(settings.uncertainty){
			targets[i].calculateTransmission(targets[i].get_ana_num_ratio());
			targets[i].calculateResonantScattering(energy_bins);
		}
	}
}

string Experiment::result_string(unsigned int n_setting) const {

	stringstream resss;
	resss<< HORIZONTAL_LINE << "\n";
	resss << ">>> SeAN RESULTS #" << n_setting << "\n";
	resss << "TARGET NAME\tRESONANT SCATTERING";
	if(settings.uncertainty)
		resss << " +- CS_UNCERTAINTY +- SCATTERING_UNCERTAINTY";
	resss << "\n";

	return resss.str();
}

string Experiment::uncertainty_string() const {

	unsigned int ntargets = (unsigned int) settings.targetNames.size();	

	stringstream unss;
	
	for(unsigned int i = 0; i < ntargets; ++i){
		unss << HORIZONTAL_LINE << "\n";
		unss << ">>> TARGET #" << i+1 << "\t: NUMERICAL UNCERTAINY ESTIMATES\n";
		unss << targets[i].uncertainty_string();
	}

	return unss.str();
};

void Experiment::print_results(unsigned int n_setting){

	unsigned int ntargets = (unsigned int) settings.targetNames.size();	

	if(settings.uncertainty){
		cout << uncertainty_string();
	}
	cout << result_string(n_setting);

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
	filename << outputfile;

	ofstream ofile;
	ofile.open(filename.str(), std::ios_base::out | std::ios_base::app);

        if(!ofile.is_open()){
                cout << "Error: " << __FILE__ << ":" << __LINE__ << ": "; 
		cout << " write_results(): File '" << outputfile << "' could not be opened." << endl;
		abort();
	}
	
	unsigned int ntargets = (unsigned int) settings.targetNames.size();	

	if(settings.uncertainty)
		ofile << uncertainty_string();
	ofile << result_string(n_setting);
	ofile.close();

	for(unsigned int i = 0; i < ntargets; ++i){
		targets[i].write_results(outputfile);
	}
}
