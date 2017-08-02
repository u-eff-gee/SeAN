#include <iostream>
#include <fstream>
#include <string>
#include <regex>

#include "InputFileReader.h"

using std::cout;
using std::endl;
using std::ifstream;
using std::string;
using std::regex;
using std::regex_replace;

void InputFileReader::readInputFile(const char* filename){
	ifname = filename;

	ifstream ifile(filename);

        if(!ifile.is_open())
                cout << "Error: ReadInputFile.cc: readInputFile(): File '" << filename << "' not found." << endl;

        cout << "> Reading input file '" << filename << "'" << endl;

        string line;
	int n = 0;
	int nline = 0;
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

			start = stop + 1;

			do{
				stop = line.substr(start, line.length()).find(",");
				if(stop == string::npos)
					stop = line.length();
				beamParams.push_back(atof(line.substr(start, stop).c_str()));
				start = stop + 1;
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
				start = stop + 1;
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
				start = stop + 1;
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
				start = stop + 1;
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
				start = stop + 1;
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

			start = stop + 1;

			do{
				stop = line.substr(start, line.length()).find(",");
				if(stop == string::npos)
					stop = line.length();
				targets[n]->addVDistParameter(atof(line.substr(start, stop).c_str()));
				start = stop + 1;
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
}

void InputFileReader::print(){
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

