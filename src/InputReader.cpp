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
using std::getline;
using std::string;
using std::regex;
using std::regex_replace;

void InputReader::readFile(vector<Settings> &settings){
	
	ifstream ifile(settings[0].inputfile);

        if(!ifile.is_open()){
                cout << "Error: " << __FILE__ << ":" << __LINE__ << ": "; 
		cout << " readFile(): File '" << settings[0].inputfile << "' not found." << endl;
		abort();
	}
        cout << "> Reading input file '" << settings[0].inputfile << "'" << endl;

        string line;
	unsigned int ntarget = 0;
	unsigned int nline = 0;

	bool emin_emax_found = false;

	// Read EMIN, EMAX
	while(!emin_emax_found){

		getline(ifile, line);
		if(line.substr(0,1) == COMMENT)
			continue;

		istringstream stream(line);
		readEminEmax(stream, settings);
		emin_emax_found = true;
	}

	bool incident_beam_found = false;

	// Read incident beam
//	while(!incident_beam_found){
//
//		getline(ifile, line);
//
//		if(line.substr(0,1) == COMMENT)
//			continue;
//
//		istringstream stream(line);
//		getline(stream, value, DELIMITER);
//		if(value == "const"){
//			settings.incidentBeam= incidentBeamModel::constant;
//			getline(stream, value, DELIMITER);
//			settings.incidentBeamParams.push_back(atof(value.c_str()));
//			incident_beam_found = true;
//		}
//
//		else if(value == "gauss"){
//			settings.incidentBeam= incidentBeamModel::gauss;
//			getline(stream, value, DELIMITER);
//			settings.incidentBeamParams.push_back(atof(value.c_str()));
//			getline(stream, value, DELIMITER);
//			settings.incidentBeamParams.push_back(atof(value.c_str()));
//			getline(stream, value, DELIMITER);
//			settings.incidentBeamParams.push_back(atof(value.c_str()));
//			incident_beam_found = true;
//		} 
//
//		else if(value == "arb"){
//			settings.incidentBeam= incidentBeamModel::arb;
//			getline(stream, value, DELIMITER);
//			settings.incidentBeamFile = regex_replace(value, regex("\\s+"), "");				
//			incident_beam_found = true;
//		}
//		
//		else{
//			cout << "Error: " << __FILE__ << ":" << __LINE__ << ": "; 
//			cout << " readFile(): Unknown option '" << value << "' for incident beam." << endl;
//			abort();
//		}
//
//	}
//
//	bool nbins_e_found = false;
//
//	// Read NBINS_ENERGY
//	while(!nbins_e_found){
//
//		getline(ifile, line);
//
//		if(line.substr(0,1) == COMMENT)
//			continue;
//
//		istringstream stream(line);
//		getline(stream, value, DELIMITER);
//		settings.nbins_e = (unsigned int) atoi(value.c_str());
//
//		nbins_e_found = true;
//	}
//
//	bool nbins_z_found = false;
//
//	// Read NBINS_Z
//	while(!nbins_z_found){
//
//		getline(ifile, line);
//
//		if(line.substr(0,1) == COMMENT)
//			continue;
//
//		istringstream stream(line);
//		getline(stream, value, DELIMITER);
//		settings.nbins_z = (unsigned int) atoi(value.c_str());
//
//		nbins_z_found = true;
//	}
//
//	// Read information about targets
//        while(getline(ifile, line)){
//		if(line.substr(0,1) == COMMENT)
//			continue;
//
//		istringstream stream(line);
//
//		switch(nline % N_TARGET_SETTINGS){
//			// Read target name
//			case 0:
//				settings.targetNames.push_back(stream.str());
//				++nline;
//				break;
//			
//			// Read resonance energies
//			case 1:
//				settings.energy.push_back(vector<double>());
//				while(getline(stream, value, DELIMITER)){
//					settings.energy[ntarget].push_back(atof(value.c_str()));
//				}
//				++nline;
//				break;
//			// Read Gamma0
//			case 2:
//				settings.gamma0.push_back(vector<double>());
//				while(getline(stream, value, DELIMITER)){
//					settings.gamma0[ntarget].push_back(atof(value.c_str()));
//				}
//				++nline;
//				break;
//
//			// Read Gamma
//			case 3:
//				settings.gamma.push_back(vector<double>());
//				while(getline(stream, value, DELIMITER)){
//					settings.gamma[ntarget].push_back(atof(value.c_str()));
//				}
//				++nline;
//				break;
//
//			// Read Ji
//			case 4:
//				while(getline(stream, value, DELIMITER)){
//					settings.ji.push_back(atof(value.c_str()));
//				}
//				++nline;
//				break;
//
//			// Read Jj
//			case 5:
//				settings.jj.push_back(vector<double>());
//				while(getline(stream, value, DELIMITER)){
//					settings.jj[ntarget].push_back(atof(value.c_str()));
//				}
//				++nline;
//				break;
//
//			// Read velocity distribution
//			case 6:
//				settings.vDistParams.push_back(vector<double>());
//				getline(stream, value, DELIMITER);
//				if(value == "zero"){
//					settings.vDist.push_back(vDistModel::zero);
//					settings.vDistFile.push_back("");
//					getline(stream, value, DELIMITER);
//					settings.vDistParams[ntarget].push_back(0.);
//				}
//				else if(value == "arb"){
//					settings.vDist.push_back(vDistModel::arb);
//					getline(stream, value, DELIMITER);
//					settings.vDistFile.push_back(regex_replace(value, regex("\\s+"), ""));
//				}
//				else if(value == "maxwell_boltzmann"){
//					settings.vDist.push_back(vDistModel::mb);
//					settings.vDistFile.push_back("");
//					getline(stream, value, DELIMITER);
//					settings.vDistParams[ntarget].push_back(atof(value.c_str()));
//				}
//				else if(value == "maxwell_boltzmann_approximation"){
//					settings.vDist.push_back(vDistModel::mba);
//					settings.vDistFile.push_back("");
//					getline(stream, value, DELIMITER);
//					settings.vDistParams[ntarget].push_back(atof(value.c_str()));
//				}
//				else if(value == "maxwell_boltzmann_debye"){
//					settings.vDist.push_back(vDistModel::mbd);
//					settings.vDistFile.push_back("");
//					getline(stream, value, DELIMITER);
//					settings.vDistParams[ntarget].push_back(atof(value.c_str()));
//					getline(stream, value, DELIMITER);
//					settings.vDistParams[ntarget].push_back(atof(value.c_str()));
//				}
//				else if(value == "maxwell_boltzmann_approximation_debye"){
//					settings.vDist.push_back(vDistModel::mbad);
//					settings.vDistFile.push_back("");
//					getline(stream, value, DELIMITER);
//					settings.vDistParams[ntarget].push_back(atof(value.c_str()));
//					getline(stream, value, DELIMITER);
//					settings.vDistParams[ntarget].push_back(atof(value.c_str()));
//				}
//
//				else{
//					cout << "Error: " << __FILE__ << ":" << __LINE__ << ": "; 
//					cout << " readFile(): Unknown option '" << value << "' for velocity distribution." << endl;
//					abort();
//				}
//				++nline;
//				break;
//
//			// Read atomic mass
//			case 7:
//				if(regex_search(line, regex("[a-zA-Z]"))){
//					settings.mass.push_back(readAME(line));
//				} else{
//					settings.mass.push_back(atof(line.c_str()));
//				}
//				++nline;
//				break;
//
//			// Read mass attenuation
//			case 8:
//				settings.mAttParams.push_back(vector<double>());
//				getline(stream, value, DELIMITER);
//				if(value == "const"){
//					settings.mAtt.push_back(mAttModel::constant);
//					getline(stream, value, DELIMITER);
//					settings.mAttParams[ntarget].push_back(atof(value.c_str()));
//					settings.mAttFile.push_back("");
//				}
//				else if(value == "nist"){
//					settings.mAtt.push_back(mAttModel::nist);
//					getline(stream, value, DELIMITER);
//					settings.mAttParams[ntarget].push_back(0.);
//					settings.mAttFile.push_back(regex_replace(value, regex("\\s+"), ""));
//				}
//				else if(value == "arb"){
//					settings.mAtt.push_back(mAttModel::arb);
//					getline(stream, value, DELIMITER);
//					settings.mAttParams[ntarget].push_back(0.);
//					settings.mAttFile.push_back(regex_replace(value, regex("\\s+"), ""));
//				}
//				else{
//					cout << "Error: " << __FILE__ << ":" << __LINE__ << ": "; 
//					cout << " readFile(): Unknown option '" << value << "' for mass attenuation." << endl;
//					abort();
//				}
//				++nline;
//				break;
//
//			// Read target thickness
//			case 9:
//				getline(stream, value, DELIMITER);
//				settings.thickness.push_back(atof(value.c_str()));
//				++nline;
//				break;
//
//			// Read target velocity
//			case 10:
//				getline(stream, value, DELIMITER);
//				settings.velocity.push_back(atof(value.c_str()));
//				++nline;
//				++ntarget;
//				break;
//		}
//        }

	ifile.close();
}

int InputReader::readDoubles(vector<double> &values, const string &value_string){
	long unsigned int loop_start = value_string.find(LOOP_START);	
	long unsigned int sim_loop_start = value_string.find(SIM_LOOP_START);	

	double start = 0.;
       	double stop = 0.;
	double inc = 0.;
	int n = 1;
	long unsigned int del1 = 0;
	long unsigned int del2 = 0;
	long unsigned int length = value_string.length();

	if(loop_start != string::npos){

		del1 = value_string.find(LOOP_DELIMITER);
		start = atof(value_string.substr(loop_start + 1, del1 - loop_start).c_str());
		
		del2 = value_string.substr(del1 + 1, length - del1).find(LOOP_DELIMITER);
		stop = atof(value_string.substr(del1 + 1, del2).c_str());

		del1 = del1 + del2 + 1;
		del2 = value_string.substr(del1,  length - del1).find(LOOP_STOP);
		n = atoi(value_string.substr(del1 + 1, del2).c_str());

		inc = (stop - start)/(n - 1);

		for(int i = 0; i < n; ++i){
			values.push_back(start + i*inc);
		}
		return 1;
	}

	else if(sim_loop_start != string::npos){

		del1 = value_string.find(LOOP_DELIMITER);
		start = atof(value_string.substr(loop_start + 1, del1 - loop_start).c_str());

		del2 = value_string.substr(del1 + 1, length - del1).find(LOOP_DELIMITER);
		stop = atof(value_string.substr(del1 + 1, del2).c_str());

		del1 = del1 + del2 + 1;
		del2 = value_string.substr(del1,  length - del1).find(SIM_LOOP_STOP);
		n = atoi(value_string.substr(del1 + 1, del2).c_str());

		inc = (stop - start)/(n - 1);

		cout << start << ", " << stop << ", " << n << endl;

		for(int i = 0; i < n; ++i){
			values.push_back(start + i*inc);
		}
		return 2;
	} else{
		values.push_back(atof(value_string.c_str()));
		return 0;
	}
}

void InputReader::readEminEmax(istringstream &stream, vector<Settings> &settings){
	vector<double> values;
	string value_string = "";
	int flag = 0;
	long unsigned int n_settings = 0;
	long unsigned int n_values = 0;

	// Read emin
	getline(stream, value_string, DELIMITER);
	flag = readDoubles(values, value_string);

	if(flag == 0){
		for(unsigned int i = 0; i < values.size(); ++i){
			settings[i].emin = values[i];
		}
	}
	if(flag == 1){
		n_settings = settings.size();
		for(unsigned int i = 0; i < values.size(); ++i){
			if(i > 0){
				for(unsigned int j = 0; j < n_settings; ++j){
					settings.push_back(settings[j]);
				}
			}
			for(unsigned int j = 0; j < n_settings; ++j){
				settings[i*n_settings + j].emin = values[i];
			}
		}
	}
	if(flag == 2){
		n_settings = settings.size();
		n_values = values.size();

		if(n_settings == 1){
			for(unsigned int i = 0; i < n_values - 1; ++i){
				settings.push_back(settings[0]);
			}
		}

		for(unsigned int i = 0; i < n_values; ++i){
			settings[i].emin = values[i];
		}
	}

	values.clear();

	// Read emax
	getline(stream, value_string, DELIMITER);
	flag = readDoubles(values, value_string);

	if(flag == 0){
		for(unsigned int i = 0; i < settings.size(); ++i){
			settings[i].emax = values[0];
		}
	}
	if(flag == 1){
		n_settings = settings.size();
		for(unsigned int i = 0; i < values.size(); ++i){
			if(i > 0){
				for(unsigned int j = 0; j < n_settings; ++j){
					settings.push_back(settings[j]);
				}
			}
			for(unsigned int j = 0; j < n_settings; ++j){
				settings[i*n_settings + j].emax = values[i];
			}
		}
	}
	if(flag == 2){
		n_settings = settings.size();
		n_values = values.size();

		if(n_settings == 1){
			for(unsigned int i = 0; i < n_values - 1; ++i){
				settings.push_back(settings[0]);
			}
		}

		for(unsigned int i = 0; i < n_values; ++i){
			settings[i].emax = values[i];
		}
	}
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
	filename << MASS_DIR << AME_FILE_NAME;
	ifstream ifile(filename.str());	

        if(!ifile.is_open()){
		cout << "Error: " << __FILE__ << ":" << __LINE__ << ": "; 
		cout << " readAME(): File '" << filename.str() << "' not found." << endl;
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

void InputReader::readNIST(vector< vector<double> > &mass_attenuation_file, const string mass_attenuation_filename){

	ifstream ifile(mass_attenuation_filename);

        if(!ifile.is_open()){
		cout << "Error: " << __FILE__ << ":" << __LINE__ << ": "; 
		cout << " readNIST(): File '" << mass_attenuation_filename << "' not found." << endl;
		abort();
	}
        cout << "> Reading input file '" << mass_attenuation_filename << "'" << endl;

	string line;

	while(getline(ifile, line)){
		if(line.substr(0,1) == COMMENT)
			continue;

		// Ignore lines with x-ray resonances since they have the same energy value as the previous bins. Those steps can not be interpolated.
		if(regex_replace(line.substr(NIST_XRAY, NIST_XRAY_LENGTH), regex("\\s+"), "") != "")
			continue;	
		
		mass_attenuation_file[0].push_back(atof(line.substr(NIST_ENERGY, NIST_ENERGY_LENGTH).c_str()));
		mass_attenuation_file[1].push_back(atof(line.substr(NIST_MU, NIST_MU_LENGTH).c_str()));
	}
}

void InputReader::read2ColumnFile(vector< vector<double> > &data, string filename){

	string line, value;
	ifstream ifile;

	ifile.open(filename);	

        if(!ifile.is_open()){
		cout << "Error: " << __FILE__ << ":" << __LINE__ << ": "; 
		cout << " read2ColumnFile(): File '" << filename << "' not found." << endl;
		abort();
	}
        cout << "> Reading input file '" << filename << "'" << endl;

	data.push_back(vector<double>());
	data.push_back(vector<double>());
	
	while(getline(ifile, line)){
		if(line.substr(0,1) == COMMENT)
			continue;

		istringstream stream(line);

		getline(stream, value, DELIMITER);
		data[0].push_back(atof(value.c_str()));
		getline(stream, value, DELIMITER); 
		data[1].push_back(atof(value.c_str()));
	}
}
