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

void Experiment::readInputFile(const char* filename){
	ifname = filename;

	ifstream ifile(filename);

        if(!ifile.is_open()){
                cout << "Error: Experiment.cc: readInputFile(): File '" << filename << "' not found." << endl;
		abort();
	}
        cout << "> Reading input file '" << filename << "'" << endl;

        string line;
	unsigned int n = 0;
	unsigned int nline = 0;
	size_t start, stop;

        while(getline(ifile, line)){
		if(line.substr(0,1) == COMMENT)
			continue;

		if(nline == 0){
			start = 0;
			stop = line.find(","); 
			emin = atof(line.substr(start, stop).c_str());
			start = stop + 1;
			stop = line.length() - 1;
			emax = atof(line.substr(start, stop).c_str());

			++nline;
			continue;
		}

		if(nline == 1){
			start = 0;
			stop = line.substr(start, line.length() - 1).find(",");
			if(stop == string::npos)
				stop = line.length();

			beam_ID = regex_replace(line.substr(start, stop), regex("\\s+"), "");

			start += stop + 1;

			do{
				stop = line.substr(start, line.length()).find(",");
				if(stop == string::npos)
					stop = line.length();
				beamParams.push_back(atof(line.substr(start, stop).c_str()));
				start += stop + 1;
			} while( stop != line.length());
			++nline;
			continue;
		}

		if(nline - INPUT_HEADER == n*NPARAMETERS + NAME){
			targets.push_back(new Target(line, n));
			++nline;
			continue;
		}

		if(nline - INPUT_HEADER == n*NPARAMETERS + E0){
			start = 0;
			do{
				stop = line.substr(start, line.length()).find(",");
				if(stop == string::npos)
					stop = line.length();
				targets[n]->addEnergy(atof(line.substr(start, stop).c_str()));
				start += stop + 1;
			} while( stop != line.length());
			++nline;
			continue;
		}

		if(nline - INPUT_HEADER == n*NPARAMETERS + GAMMA0){
			start = 0;
			do{
				stop = line.substr(start, line.length()).find(",");
				if(stop == string::npos)
					stop = line.length();
				targets[n]->addGamma0(atof(line.substr(start, stop).c_str()));
				start += stop + 1;
			} while( stop != line.length());
			++nline;
			continue;
		}

		if(nline - INPUT_HEADER == n*NPARAMETERS + GAMMA){
			start = 0;
			do{
				stop = line.substr(start, line.length()).find(",");
				if(stop == string::npos)
					stop = line.length();
				targets[n]->addGamma(atof(line.substr(start, stop).c_str()));
				start += stop + 1;
			} while( stop != line.length());
			++nline;
			continue;
		}

		if(nline - INPUT_HEADER == n*NPARAMETERS + JJ){
			start = 0;
			do{
				stop = line.substr(start, line.length()).find(",");
				if(stop == string::npos)
					stop = line.length();
				targets[n]->addJJ(atof(line.substr(start, stop).c_str()));
				start += stop + 1;
			} while( stop != line.length());
			++nline;
			continue;
		}

		if(nline - INPUT_HEADER == n*NPARAMETERS + J0){
			targets[n]->setJ0(atof(line.c_str()));
			++nline;
			continue;
		}
	
		if(nline - INPUT_HEADER == n*NPARAMETERS + V){
			start = 0;
			stop = line.substr(start, line.length() - 1).find(",");
			if(stop == string::npos)
				stop = line.length();

			targets[n]->setVDist(regex_replace(line.substr(start, stop), regex("\\s+"), ""));

			start += stop + 1;

			do{
				stop = line.substr(start, line.length()).find(",");
				if(stop == string::npos)
					stop = line.length();
				targets[n]->addVDistParameter(atof(line.substr(start, stop).c_str()));
				start += stop + 1;
			} while( stop != line.length());
			++nline;
			continue;
		}
		
		if(nline - INPUT_HEADER == n*NPARAMETERS + M){
			targets[n]->setMass(atof(line.c_str()));
			++nline;
			continue;
		}

		if(nline - INPUT_HEADER == n*NPARAMETERS + MU){
			targets[n]->setMassAttenuation(regex_replace(line, regex("\\s+"), ""));
			++nline;
			continue;
		}

		if(nline - INPUT_HEADER == n*NPARAMETERS + Z){
			targets[n]->setZ(atof(line.c_str()));
			++nline;
			++n;
			continue;
		}

        }

	ifile.close();

	createEnergyBins(emin, emax);
}

void Experiment::testIntegration(bool plot){
	createEnergyBins(emin, emax);
	targets[0]->testIntegration(emin, emax, energy_bins, beamParams);
	if(plot){
		targets[0]->plotTestIntegration(energy_bins);
	}
}

void Experiment::createEnergyBins(double emin, double emax){
	double delta_e = (emax - emin)/NBINS;

	for(int i = 0; i < NBINS; ++i)
		energy_bins[i] = emin + i*delta_e;
};

void Experiment::crossSections(bool plot){
	for(unsigned int i = 0; i < targets.size(); ++i){
		targets[i]->calculateCrossSection(energy_bins);
		targets[i]->calculateVelocityDistribution(energy_bins);
		targets[i]->calculateDopplerShift(energy_bins);

		if(plot){
			targets[i]->plotCrossSection(energy_bins);
			targets[i]->plotVelocityDistribution();
			targets[i]->plotDopplerShift(energy_bins);
		}
	}
};

void Experiment::transmission(bool plot){
	targets[0]->calculateIncidentBeam(energy_bins, beam_ID, beamParams);
	for(unsigned int i = 0; i < targets.size(); ++i){
		targets[i]->calculateMassAttenuation(energy_bins);
		targets[i]->calculateZBins();
		if(i > 0){
			targets[i]->setIncidentBeam(targets[i-1]->getTransmittedBeam());
		}
		targets[i]->calculatePhotonFluxDensity();
		if(i == 0){
			targets[i]->calculateTransmittedBeam();
		}
		targets[i]->calculateResonanceAbsorptionDensity();
		targets[i]->calculateAbsorption(energy_bins);

		if(plot){
			targets[i]->plotMassAttenuation(energy_bins);
			targets[i]->plotMu();
			targets[i]->plotPhotonFluxDensity(energy_bins);
			targets[i]->plotResonanceAbsorptionDensity(energy_bins);
		}
	}
}

void Experiment::print(){
	cout << "INPUT FILE '" << ifname << "'" << endl;
	cout << "[EMIN, EMAX] = [" << emin << ", " << emax << "]" << endl;
	cout << "BEAM = " << beam_ID;

	if( beamParams.size() ){
		cout << " ( ";
		for(unsigned int i = 0; i < beamParams.size(); ++i)
			cout << beamParams[i] << " ";

		cout << ")" << endl;
	} else{
		cout << endl;
	}

	for(unsigned int i = 0; i < targets.size(); ++i){
		cout << endl;
		targets[i]->print();
	}
}

