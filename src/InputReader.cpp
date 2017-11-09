#include "InputReader.h"
#include "Config.h"

#include <iostream>
#include <fstream>
#include <string>
#include <regex>
#include <sstream>

using std::cout;
using std::endl;
using std::ifstream;
using std::stringstream;
using std::istringstream;
using std::getline;
using std::string;
using std::regex;
using std::regex_replace;

void InputReader::readFile(Settings &settings){
	
	ifstream ifile(settings.inputfile);

        if(!ifile.is_open()){
                cout << "Error: " << __FILE__ << ":" << __LINE__ << ": "; 
		cout << " readFile(): File '" << settings.inputfile << "' not found." << endl;
		abort();
	}
        cout << "> Reading input file '" << settings.inputfile << "'" << endl;

        string line, value;
	unsigned int ntarget = 0;
	unsigned int nline = 0;
	//size_t start, stop;

        while(getline(ifile, line)){
		if(line.substr(0,1) == COMMENT)
			continue;

		istringstream stream(line);

		switch(nline % (N_EXPERIMENT_SETTINGS + N_TARGET_SETTINGS)){
			case 0:
				getline(stream, value, DELIMITER);
				settings.emin = atof(value.c_str());
				getline(stream, value, DELIMITER); 
				settings.emax = atof(value.c_str());
				++nline;
				break;

			case 1:
				getline(stream, value, DELIMITER);
				if(value == "const"){
					settings.incidentBeam= incidentBeamModel::constant;
					getline(stream, value, DELIMITER);
					settings.incidentBeamParams.push_back(atof(value.c_str()));
				}

				else if(value == "gauss"){
					settings.incidentBeam= incidentBeamModel::gauss;
					getline(stream, value, DELIMITER);
					settings.incidentBeamParams.push_back(atof(value.c_str()));
					getline(stream, value, DELIMITER);
					settings.incidentBeamParams.push_back(atof(value.c_str()));
					getline(stream, value, DELIMITER);
					settings.incidentBeamParams.push_back(atof(value.c_str()));
				} 

				else if(value == "arb"){
					settings.incidentBeam= incidentBeamModel::arb;
					getline(stream, value, DELIMITER);
					settings.incidentBeamFile = value;				
				}
				
				else{
					cout << "Error: " << __FILE__ << ":" << __LINE__ << ": "; 
					cout << " readFile(): Unknown option '" << value << "' for incident beam." << endl;
					abort();
				}
				++nline;
				break;

			case 2:
				getline(stream, value, DELIMITER);
				settings.nbins_e = (unsigned int) atoi(value.c_str());
				++nline;
				break;

			case 3:
				getline(stream, value, DELIMITER);
				settings.nbins_z = (unsigned int) atoi(value.c_str());
				++nline;
				break;

			case 4:
				settings.targetNames.push_back(stream.str());
				++nline;
				break;
			
			case 5:
				settings.energy.push_back(vector<double>());
				while(getline(stream, value, DELIMITER)){
					settings.energy[ntarget].push_back(atof(value.c_str()));
				}
				++nline;
				break;

			case 6:
				settings.gamma0.push_back(vector<double>());
				while(getline(stream, value, DELIMITER)){
					settings.gamma0[ntarget].push_back(atof(value.c_str()));
				}
				++nline;
				break;

			case 7:
				settings.gamma.push_back(vector<double>());
				while(getline(stream, value, DELIMITER)){
					settings.gamma[ntarget].push_back(atof(value.c_str()));
				}
				++nline;
				break;

			case 8:
				settings.ji.push_back(vector<double>());
				while(getline(stream, value, DELIMITER)){
					settings.ji[ntarget].push_back(atof(value.c_str()));
				}
				++nline;
				break;

			case 9:
				settings.jj.push_back(vector<double>());
				while(getline(stream, value, DELIMITER)){
					settings.jj[ntarget].push_back(atof(value.c_str()));
				}
				++nline;
				break;

			case 10:
				getline(stream, value, DELIMITER);
				if(value == "zero"){
					settings.vDist.push_back(vDistModel::zero);
					getline(stream, value, DELIMITER);
					settings.incidentBeamParams.push_back(atof(value.c_str()));
				}
				else if(value == "maxwell_boltzmann"){
					settings.vDist.push_back(vDistModel::mb);
					getline(stream, value, DELIMITER);
					settings.incidentBeamParams.push_back(atof(value.c_str()));
				}
				else if(value == "maxwell_boltzmann_approximation"){
					settings.vDist.push_back(vDistModel::mba);
					getline(stream, value, DELIMITER);
					settings.incidentBeamParams.push_back(atof(value.c_str()));
				}

				else{
					cout << "Error: " << __FILE__ << ":" << __LINE__ << ": "; 
					cout << " readFile(): Unknown option '" << value << "' for velocity distribution." << endl;
					abort();
				}
				++nline;
				break;

			case 11:
				if(regex_search(line, regex("[a-zA-Z]"))){
					settings.mass.push_back(readAME(line));
				} else{
					settings.mass.push_back(atof(line.c_str()));
				}
				++nline;
				break;

			case 12:
				getline(stream, value, DELIMITER);
				if(value == "const"){
					settings.mAtt.push_back(mAttModel::constant);
					getline(stream, value, DELIMITER);
					settings.mAttParams.push_back(atof(value.c_str()));
				}
				else if(value == "nist"){
					settings.mAtt.push_back(mAttModel::nist);
					getline(stream, value, DELIMITER);
					settings.mAttFile = value;
				}
				else if(value == "arb"){
					settings.mAtt.push_back(mAttModel::arb);
					getline(stream, value, DELIMITER);
					settings.mAttFile = value;
				}
				else{
					cout << "Error: " << __FILE__ << ":" << __LINE__ << ": "; 
					cout << " readFile(): Unknown option '" << value << "' for mass attenuation." << endl;
					abort();
				}
				++nline;
				break;

			case 13:
				getline(stream, value, DELIMITER);
				settings.thickness.push_back(atof(value.c_str()));
				++nline;
				break;

			case 14:
				getline(stream, value, DELIMITER);
				settings.velocity.push_back(atof(value.c_str()));
				++nline;
				break;
		}
        }

	ifile.close();
}

double InputReader::readAME(string isotope){
	unsigned int separator = 0;

	for(unsigned int i = 1; i <= isotope.length(); ++i){
		if(regex_search(isotope.substr(0,i), regex("[a-zA-Z]"))){
			separator = i - 1;
			break;
		}
	}

	int mass_number = atoi(isotope.substr(0,separator).c_str());
	string string_mass_number = isotope.substr(0,separator);
	string isotope_name = isotope.substr(separator, isotope.length() - separator + 1);

	stringstream filename;
	filename << MASS_DIR << "mass_list.txt";
	ifstream ifile(filename.str());	

        if(!ifile.is_open()){
                cout << "Error: Target.cc: readAME(): File '" << filename.str() << "' not found." << endl;
		abort();
	}
        cout << "> Reading input file '" << filename.str() << "'" << endl;

        string line;
	unsigned int nline = 0;

        while(getline(ifile, line)){
		if(nline > AME_HEADER_LENGTH){
			if(atoi(line.substr(AME_MASS_NUMBER, 3).c_str()) == mass_number && regex_replace(line.substr(AME_ISOTOPE, 2), regex("\\s+"), "") == isotope_name){
				return atof(regex_replace(line.substr(AME_MASS_START, AME_MASS_LENGTH), regex("\\s+"), "").c_str())*1.0e-6;
			}
		}
		++nline;
	}

	return 0.;
}

