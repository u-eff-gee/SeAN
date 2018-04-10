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
        //cout << "> Reading input file '" << settings[0].inputfile << "'" << endl;

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
	while(!incident_beam_found){

		getline(ifile, line);

		if(line.substr(0,1) == COMMENT)
			continue;

		istringstream stream(line);
		readIncidentBeam(stream, settings);

		incident_beam_found = true;
	}

	bool nbins_e_found = false;

	// Read NBINS_E
	while(!nbins_e_found){

		getline(ifile, line);

		if(line.substr(0,1) == COMMENT)
			continue;

		istringstream stream(line);
		readNBins_E(stream, settings);

		nbins_e_found = true;
	}

	bool nbins_z_found = false;

	// Read NBINS_Z
	while(!nbins_z_found){

		getline(ifile, line);

		if(line.substr(0,1) == COMMENT)
			continue;

		istringstream stream(line);
		readNBins_Z(stream, settings);

		nbins_z_found = true;
	}

	// Read information about targets
        while(getline(ifile, line)){
		if(line.substr(0,1) == COMMENT)
			continue;

		istringstream stream(line);

		switch(nline % N_TARGET_SETTINGS){
			// Read target name
			case 0:
				for(unsigned i = 0; i < settings.size(); ++i){
					settings[i].targetNames.push_back(stream.str());
				}
				++nline;
				break;
			
			// Read resonance energies
			case 1:
				for(unsigned i = 0; i < settings.size(); ++i){
					settings[i].energy.push_back(vector<double>());
				}
				readEnergy(stream, settings, ntarget);
				++nline;
				break;
			// Read Gamma0
			case 2:
				for(unsigned i = 0; i < settings.size(); ++i){
					settings[i].gamma0.push_back(vector<double>());
				}
				readGamma0(stream, settings, ntarget);
				++nline;
				break;

			// Read Gamma
			case 3:
				for(unsigned i = 0; i < settings.size(); ++i){
					settings[i].gamma.push_back(vector<double>());
				}
				readGamma(stream, settings, ntarget);
				++nline;
				break;

			// Read J0
			case 4:
				readJ0(stream, settings, ntarget);
				++nline;
				break;
				
			// Read Jj
			case 5:
				for(unsigned i = 0; i < settings.size(); ++i){
					settings[i].jj.push_back(vector<double>());
				}
				readJ(stream, settings, ntarget);
				++nline;
				break;
				
			// Read doppler broadening
			case 6:
				for(unsigned i = 0; i < settings.size(); ++i){
					settings[i].dopplerParams.push_back(vector<double>());
				}
				readDopplerBroadening(stream, settings, ntarget);
				++nline;
				break;

			// Read atomic mass
			case 7:
				readMass(stream, settings);
				++nline;
				break;

			// Read mass attenuation
			case 8:
				for(unsigned i = 0; i < settings.size(); ++i){
					settings[i].mAttParams.push_back(vector<double>());
				}
				readMassAttenuation(stream, settings, ntarget);
				++nline;
				break;

			// Read target thickness
			case 9:
				readTargetThickness(stream, settings);
				++nline;
				break;
				
			// Read target velocity
			case 10:
				readVelocity(stream, settings);
				++nline;
				++ntarget;
				break;
		}
        }

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
		start = atof(value_string.substr(sim_loop_start + 1, del1 - loop_start).c_str());

		del2 = value_string.substr(del1 + 1, length - del1).find(LOOP_DELIMITER);
		stop = atof(value_string.substr(del1 + 1, del2).c_str());

		del1 = del1 + del2 + 1;
		del2 = value_string.substr(del1,  length - del1).find(SIM_LOOP_STOP);
		n = atoi(value_string.substr(del1 + 1, del2).c_str());

		inc = (stop - start)/(n - 1);

		for(int i = 0; i < n; ++i){
			values.push_back(start + i*inc);
		}
		return 2;
	} else{
		values.push_back(atof(value_string.c_str()));
		return 0;
	}
}

int InputReader::readInts(vector<int> &values, const string &value_string){
	long unsigned int loop_start = value_string.find(LOOP_START);	
	long unsigned int sim_loop_start = value_string.find(SIM_LOOP_START);	

	int start = 0;
       	int stop = 0;
	int inc = 0;
	int n = 1;
	long unsigned int del1 = 0;
	long unsigned int del2 = 0;
	long unsigned int length = value_string.length();

	if(loop_start != string::npos){

		del1 = value_string.find(LOOP_DELIMITER);
		start = atoi(value_string.substr(loop_start + 1, del1 - loop_start).c_str());
		
		del2 = value_string.substr(del1 + 1, length - del1).find(LOOP_DELIMITER);
		stop = atoi(value_string.substr(del1 + 1, del2).c_str());

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
		start = atoi(value_string.substr(sim_loop_start + 1, del1 - loop_start).c_str());

		del2 = value_string.substr(del1 + 1, length - del1).find(LOOP_DELIMITER);
		stop = atoi(value_string.substr(del1 + 1, del2).c_str());

		del1 = del1 + del2 + 1;
		del2 = value_string.substr(del1,  length - del1).find(SIM_LOOP_STOP);
		n = atoi(value_string.substr(del1 + 1, del2).c_str());

		inc = (stop - start)/(n - 1);

		for(int i = 0; i < n; ++i){
			values.push_back(start + i*inc);
		}
		return 2;
	} else{
		values.push_back(atoi(value_string.c_str()));
		return 0;
	}
}

void InputReader::readEminEmax(istringstream &stream, vector<Settings> &settings){
	vector<double> values;
	string value_string = "";
	int flag = 0;
	long unsigned int n_settings = settings.size();
	long unsigned int n_values = 0;

	// Read emin
	getline(stream, value_string, DELIMITER);
	flag = readDoubles(values, value_string);

	if(flag == 0){
		n_settings = settings.size();

		for(unsigned int i = 0; i < n_settings; ++i){
			settings[i].emin = values[0];
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
		n_values = values.size();
		n_settings = settings.size();

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
		n_settings = settings.size();
		for(unsigned int i = 0; i < n_settings; ++i){
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

void InputReader::readIncidentBeam(istringstream &stream, vector<Settings> &settings){

	vector<double> values;
	string value_string = "";
	int flag = 0;
	long unsigned int n_settings = settings.size();
	long unsigned int n_values = 0;

	// Read incident beam model
	getline(stream, value_string, DELIMITER);

	if(value_string=="const"){
		for(unsigned int i = 0; i < n_settings; ++i){
			settings[i].incidentBeam = incidentBeamModel::constant;
		}
	}
	else if(value_string=="gauss"){
		for(unsigned int i = 0; i < n_settings; ++i){
			settings[i].incidentBeam = incidentBeamModel::gauss;
		}
	}
	else if(value_string=="arb"){
		for(unsigned int i = 0; i < n_settings; ++i){
			settings[i].incidentBeam = incidentBeamModel::arb;
		}
	} else{
		cout << "Error: " << __FILE__ << ":" << __LINE__ << ": "; 
		cout << " readFile(): Unknown option '" << value_string << "' for incident beam." << endl;
		abort();
	}

	// Read parameters

	while(getline(stream, value_string, DELIMITER)){
		flag = readDoubles(values, value_string);

		if(flag == 0){
			n_settings = settings.size();
			for(unsigned int i = 0; i < n_settings; ++i){
				settings[i].incidentBeamParams.push_back(values[0]);
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
			}
			for(unsigned int i = 0; i < values.size(); ++i){
				for(unsigned int j = 0; j < n_settings; ++j){
					settings[i*n_settings + j].incidentBeamParams.push_back(values[i]);
				}
			}
		}
		if(flag == 2){
			n_values = values.size();

			if(n_settings == 1){
				for(unsigned int i = 0; i < n_values - 1; ++i){
					settings.push_back(settings[0]);
				}
			}

			for(unsigned int i = 0; i < n_values; ++i){
				settings[i].incidentBeamParams.push_back(values[i]);
			}
		}

		values.clear();
	}
}

void InputReader::readNBins_E(istringstream &stream, vector<Settings> &settings){
	vector<int> values;
	string value_string = "";
	int flag = 0;
	long unsigned int n_settings = settings.size();
	long unsigned int n_values = 0;

	// Read nbins_e
	getline(stream, value_string, DELIMITER);
	flag = readInts(values, value_string);

	if(flag == 0){
		n_settings = settings.size();

		for(unsigned int i = 0; i < n_settings; ++i){
			settings[i].nbins_e = (unsigned int) values[0];
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
				settings[i*n_settings + j].nbins_e = (unsigned int) values[i];
			}
		}
	}
	if(flag == 2){
		n_values = values.size();
		n_settings = settings.size();

		if(n_settings == 1){
			for(unsigned int i = 0; i < n_values - 1; ++i){
				settings.push_back(settings[0]);
			}
		}

		for(unsigned int i = 0; i < n_values; ++i){
			settings[i].nbins_e = (unsigned int) values[i];
		}
	}
}

void InputReader::readNBins_Z(istringstream &stream, vector<Settings> &settings){
	vector<int> values;
	string value_string = "";
	int flag = 0;
	long unsigned int n_settings = settings.size();
	long unsigned int n_values = 0;

	// Read nbins_z
	getline(stream, value_string, DELIMITER);
	flag = readInts(values, value_string);

	if(flag == 0){
		n_settings = settings.size();

		for(unsigned int i = 0; i < n_settings; ++i){
			settings[i].nbins_z = (unsigned int) values[0];
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
				settings[i*n_settings + j].nbins_z = (unsigned int) values[i];
			}
		}
	}
	if(flag == 2){
		n_values = values.size();
		n_settings = settings.size();

		if(n_settings == 1){
			for(unsigned int i = 0; i < n_values - 1; ++i){
				settings.push_back(settings[0]);
			}
		}

		for(unsigned int i = 0; i < n_values; ++i){
			settings[i].nbins_z = (unsigned int) values[i];
		}
	}
}

void InputReader::readEnergy(istringstream &stream, vector<Settings> &settings, unsigned int ntarget){

	vector<double> values;
	string value_string = "";
	int flag = 0;
	long unsigned int n_settings = settings.size();
	long unsigned int n_values = 0;

	// Read energies

	while(getline(stream, value_string, DELIMITER)){
		flag = readDoubles(values, value_string);

		if(flag == 0){
			n_settings = settings.size();
			for(unsigned int i = 0; i < n_settings; ++i){
				settings[i].energy[ntarget].push_back(values[0]);
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
			}
			for(unsigned int i = 0; i < values.size(); ++i){
				for(unsigned int j = 0; j < n_settings; ++j){
					settings[i*n_settings + j].energy[ntarget].push_back(values[i]);
				}
			}
		}
		if(flag == 2){
			n_values = values.size();

			if(n_settings == 1){
				for(unsigned int i = 0; i < n_values - 1; ++i){
					settings.push_back(settings[0]);
				}
			}

			for(unsigned int i = 0; i < n_values; ++i){
				settings[i].energy[ntarget].push_back(values[i]);
			}
		}

		values.clear();
	}
}

void InputReader::readGamma0(istringstream &stream, vector<Settings> &settings, unsigned int ntarget){

	vector<double> values;
	string value_string = "";
	int flag = 0;
	long unsigned int n_settings = settings.size();
	long unsigned int n_values = 0;

	// Read values for gamma0

	while(getline(stream, value_string, DELIMITER)){
		flag = readDoubles(values, value_string);

		if(flag == 0){
			n_settings = settings.size();
			for(unsigned int i = 0; i < n_settings; ++i){
				settings[i].gamma0[ntarget].push_back(values[0]);
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
			}
			for(unsigned int i = 0; i < values.size(); ++i){
				for(unsigned int j = 0; j < n_settings; ++j){
					settings[i*n_settings + j].gamma0[ntarget].push_back(values[i]);
				}
			}
		}
		if(flag == 2){
			n_values = values.size();

			if(n_settings == 1){
				for(unsigned int i = 0; i < n_values - 1; ++i){
					settings.push_back(settings[0]);
				}
			}

			for(unsigned int i = 0; i < n_values; ++i){
				settings[i].gamma0[ntarget].push_back(values[i]);
			}
		}

		values.clear();
	}
}

void InputReader::readGamma(istringstream &stream, vector<Settings> &settings, unsigned int ntarget){

	vector<double> values;
	string value_string = "";
	int flag = 0;
	long unsigned int n_settings = settings.size();
	long unsigned int n_values = 0;

	// Read values for gamma

	while(getline(stream, value_string, DELIMITER)){
		flag = readDoubles(values, value_string);

		if(flag == 0){
			n_settings = settings.size();
			for(unsigned int i = 0; i < n_settings; ++i){
				settings[i].gamma[ntarget].push_back(values[0]);
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
			}
			for(unsigned int i = 0; i < values.size(); ++i){
				for(unsigned int j = 0; j < n_settings; ++j){
					settings[i*n_settings + j].gamma[ntarget].push_back(values[i]);
				}
			}
		}
		if(flag == 2){
			n_values = values.size();

			if(n_settings == 1){
				for(unsigned int i = 0; i < n_values - 1; ++i){
					settings.push_back(settings[0]);
				}
			}

			for(unsigned int i = 0; i < n_values; ++i){
				settings[i].gamma[ntarget].push_back(values[i]);
			}
		}

		values.clear();
	}
}


void InputReader::readJ0(istringstream &stream, vector<Settings> &settings, unsigned int ntarget){

	vector<double> values;
	string value_string = "";
	int flag = 0;
	long unsigned int n_settings = settings.size();
	long unsigned int n_values = 0;

	// Read values for j0 

	while(getline(stream, value_string, DELIMITER)){
		flag = readDoubles(values, value_string);

		if(flag == 0){
			n_settings = settings.size();
			for(unsigned int i = 0; i < n_settings; ++i){
				settings[i].ji.push_back(values[0]);
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
			}
			for(unsigned int i = 0; i < values.size(); ++i){
				for(unsigned int j = 0; j < n_settings; ++j){
					settings[i*n_settings + j].ji.push_back(values[i]);
				}
			}
		}
		if(flag == 2){
			n_values = values.size();

			if(n_settings == 1){
				for(unsigned int i = 0; i < n_values - 1; ++i){
					settings.push_back(settings[0]);
				}
			}

			for(unsigned int i = 0; i < n_values; ++i){
				settings[i].ji.push_back(values[i]);
			}
		}

		values.clear();
	}
}

void InputReader::readJ(istringstream &stream, vector<Settings> &settings, unsigned int ntarget){

	vector<double> values;
	string value_string = "";
	int flag = 0;
	long unsigned int n_settings = settings.size();
	long unsigned int n_values = 0;

	// Read values for jj

	while(getline(stream, value_string, DELIMITER)){
		flag = readDoubles(values, value_string);

		if(flag == 0){
			n_settings = settings.size();
			for(unsigned int i = 0; i < n_settings; ++i){
				settings[i].jj[ntarget].push_back(values[0]);
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
			}
			for(unsigned int i = 0; i < values.size(); ++i){
				for(unsigned int j = 0; j < n_settings; ++j){
					settings[i*n_settings + j].jj[ntarget].push_back(values[i]);
				}
			}
		}
		if(flag == 2){
			n_values = values.size();

			if(n_settings == 1){
				for(unsigned int i = 0; i < n_values - 1; ++i){
					settings.push_back(settings[0]);
				}
			}

			for(unsigned int i = 0; i < n_values; ++i){
				settings[i].jj[ntarget].push_back(values[i]);
			}
		}

		values.clear();
	}
}

void InputReader::readDopplerBroadening(istringstream &stream, vector<Settings> &settings, const unsigned int ntarget){

	vector<double> values;
	string value_string = "";
	int flag = 0;
	long unsigned int n_settings = settings.size();
	long unsigned int n_values = 0;
	bool has_vDist_file = false;
	bool has_cs_file = false;
	bool needs_phdos = false;

	// Read velocity distribution model
	getline(stream, value_string, DELIMITER);

	if(value_string=="zero"){
		for(unsigned int i = 0; i < n_settings; ++i){
			settings[i].dopplerBroadening.push_back(dopplerModel::zero);
		}
	} else if(value_string=="arb_velocity_distribution"){
		for(unsigned int i = 0; i < n_settings; ++i){
			settings[i].dopplerBroadening.push_back(dopplerModel::arb_vdist);
		}

		has_vDist_file = true;

	} else if(value_string=="arb_cross_section"){
		for(unsigned int i = 0; i < n_settings; ++i){
			settings[i].dopplerBroadening.push_back(dopplerModel::arb_cs);
		}

		has_cs_file = true;

	} else if(value_string=="maxwell_boltzmann"){
		for(unsigned int i = 0; i < n_settings; ++i){
			settings[i].dopplerBroadening.push_back(dopplerModel::mb);
		}
	} else if(value_string=="maxwell_boltzmann_approximation"){
		for(unsigned int i = 0; i < n_settings; ++i){
			settings[i].dopplerBroadening.push_back(dopplerModel::mba);
		}
	} else if(value_string=="maxwell_boltzmann_approximation"){
		for(unsigned int i = 0; i < n_settings; ++i){
			settings[i].dopplerBroadening.push_back(dopplerModel::mba);
		}
	} else if(value_string=="maxwell_boltzmann_debye"){
		for(unsigned int i = 0; i < n_settings; ++i){
			settings[i].dopplerBroadening.push_back(dopplerModel::mbd);
		}
	} else if(value_string=="maxwell_boltzmann_approximation_debye"){
		for(unsigned int i = 0; i < n_settings; ++i){
			settings[i].dopplerBroadening.push_back(dopplerModel::mbad);
		}
	} else if(value_string=="phdos"){
		for(unsigned int i = 0; i < n_settings; ++i){
			settings[i].dopplerBroadening.push_back(dopplerModel::phdos);
		}

		needs_phdos = true;
	} else{
		cout << "Error: " << __FILE__ << ":" << __LINE__ << ": "; 
		cout << " readFile(): Unknown option '" << value_string << "' for doppler broadening." << endl;
		abort();
	}

	// Read parameters
	
	if(has_vDist_file){
		getline(stream, value_string, DELIMITER);
			
		for(unsigned int i = 0; i < n_settings; ++i){
			settings[i].velocityBinFile.push_back(regex_replace(value_string, regex("\\s+"), ""));
		}

		getline(stream, value_string, DELIMITER);
			
		for(unsigned int i = 0; i < n_settings; ++i){
			settings[i].vDistFile.push_back(regex_replace(value_string, regex("\\s+"), ""));
		}

	} else if(has_cs_file){
		getline(stream, value_string, DELIMITER);
			
		for(unsigned int i = 0; i < n_settings; ++i){
			settings[i].energyBinFile.push_back(regex_replace(value_string, regex("\\s+"), ""));
		}

		getline(stream, value_string, DELIMITER);
			
		for(unsigned int i = 0; i < n_settings; ++i){
			settings[i].crosssectionFile.push_back(regex_replace(value_string, regex("\\s+"), ""));
		}

	} else{
		if(needs_phdos){
			getline(stream, value_string, DELIMITER);
				
			for(unsigned int i = 0; i < n_settings; ++i){
				settings[i].omegaFile.push_back(regex_replace(value_string, regex("\\s+"), ""));
			}

			getline(stream, value_string, DELIMITER);
				
			for(unsigned int i = 0; i < n_settings; ++i){
				settings[i].polarizationFile.push_back(regex_replace(value_string, regex("\\s+"), ""));
			}

			getline(stream, value_string, DELIMITER);
				
			for(unsigned int i = 0; i < n_settings; ++i){
				settings[i].momentumFile.push_back(regex_replace(value_string, regex("\\s+"), ""));
			}
		}


		while(getline(stream, value_string, DELIMITER)){
			flag = readDoubles(values, value_string);

			if(flag == 0){
				n_settings = settings.size();
				for(unsigned int i = 0; i < n_settings; ++i){
					settings[i].dopplerParams[ntarget].push_back(values[0]);
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
				}
				for(unsigned int i = 0; i < values.size(); ++i){
					for(unsigned int j = 0; j < n_settings; ++j){
						settings[i*n_settings + j].dopplerParams[ntarget].push_back(values[i]);
					}
				}
			}
			if(flag == 2){
				n_values = values.size();

				if(n_settings == 1){
					for(unsigned int i = 0; i < n_values - 1; ++i){
						settings.push_back(settings[0]);
					}
				}

				for(unsigned int i = 0; i < n_values; ++i){
					settings[i].dopplerParams[ntarget].push_back(values[i]);
				}
			}

			values.clear();
		}
	}
}

void InputReader::readMass(istringstream &stream, vector<Settings> &settings){
	vector<double> values;
	string value_string = "";
	int flag = 0;
	long unsigned int n_settings = settings.size();
	long unsigned int n_values = 0;

	// Read values for the mass
	getline(stream, value_string, DELIMITER);

	if(regex_search(value_string, regex("[a-zA-Z]"))){
		values.push_back(readAME(regex_replace(value_string, regex("\\s+"), "")));
		for(unsigned int i = 0; i < n_settings; ++i){
			settings[i].mass.push_back(values[0]);
		}
	} else{
		flag = readDoubles(values, value_string);

		if(flag == 0){
			n_settings = settings.size();
			for(unsigned int i = 0; i < n_settings; ++i){
				settings[i].mass.push_back(values[0]);
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
			}
			for(unsigned int i = 0; i < values.size(); ++i){
				for(unsigned int j = 0; j < n_settings; ++j){
					settings[i*n_settings + j].mass.push_back(values[i]);
				}
			}
		}
		if(flag == 2){
			n_values = values.size();

			if(n_settings == 1){
				for(unsigned int i = 0; i < n_values - 1; ++i){
					settings.push_back(settings[0]);
				}
			}

			for(unsigned int i = 0; i < n_values; ++i){
				settings[i].mass.push_back(values[i]);
			}
		}
	}
	
	values.clear();
}

void InputReader::readMassAttenuation(istringstream &stream, vector<Settings> &settings, const unsigned int ntarget){

	vector<double> values;
	string value_string = "";
	int flag = 0;
	long unsigned int n_settings = settings.size();
	long unsigned int n_values = 0;
	bool file = false;

	// Read velocity distribution model
	getline(stream, value_string, DELIMITER);

	if(value_string=="const"){
		for(unsigned int i = 0; i < n_settings; ++i){
			settings[i].mAtt.push_back(mAttModel::constant);
		}
	} else if(value_string=="nist"){
		for(unsigned int i = 0; i < n_settings; ++i){
			settings[i].mAtt.push_back(mAttModel::nist);
		}

		file = true;

	} else if(value_string=="arb"){
		for(unsigned int i = 0; i < n_settings; ++i){
			settings[i].mAtt.push_back(mAttModel::arb);
		}

		file = true;
	} else{
		cout << "Error: " << __FILE__ << ":" << __LINE__ << ": "; 
		cout << " readFile(): Unknown option '" << value_string << "' for mass attenuation." << endl;
		abort();
	}

	// Read parameters
	
	if(file){
		getline(stream, value_string, DELIMITER);
			
		for(unsigned int i = 0; i < n_settings; ++i){
			settings[i].mAttFile.push_back(regex_replace(value_string, regex("\\s+"), ""));
		}
	} else{

		while(getline(stream, value_string, DELIMITER)){
			flag = readDoubles(values, value_string);

			if(flag == 0){
				n_settings = settings.size();
				for(unsigned int i = 0; i < n_settings; ++i){
					settings[i].mAttParams[ntarget].push_back(values[0]);
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
				}
				for(unsigned int i = 0; i < values.size(); ++i){
					for(unsigned int j = 0; j < n_settings; ++j){
						settings[i*n_settings + j].mAttParams[ntarget].push_back(values[i]);
					}
				}
			}
			if(flag == 2){
				n_values = values.size();

				if(n_settings == 1){
					for(unsigned int i = 0; i < n_values - 1; ++i){
						settings.push_back(settings[0]);
					}
				}

				for(unsigned int i = 0; i < n_values; ++i){
					settings[i].mAttParams[ntarget].push_back(values[i]);
				}
			}

			values.clear();
		}
	}
}

void InputReader::readTargetThickness(istringstream &stream, vector<Settings> &settings){
	vector<double> values;
	string value_string = "";
	int flag = 0;
	long unsigned int n_settings = settings.size();
	long unsigned int n_values = 0;

	// Read values for the target thickness
	getline(stream, value_string, DELIMITER);
	flag = readDoubles(values, value_string);

	if(flag == 0){
		n_settings = settings.size();
		for(unsigned int i = 0; i < n_settings; ++i){
			settings[i].thickness.push_back(values[0]);
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
		}
		for(unsigned int i = 0; i < values.size(); ++i){
			for(unsigned int j = 0; j < n_settings; ++j){
				settings[i*n_settings + j].thickness.push_back(values[i]);
			}
		}
	}
	if(flag == 2){
		n_values = values.size();

		if(n_settings == 1){
			for(unsigned int i = 0; i < n_values - 1; ++i){
				settings.push_back(settings[0]);
			}
		}

		for(unsigned int i = 0; i < n_values; ++i){
			settings[i].thickness.push_back(values[i]);
		}
	}
}

void InputReader::readVelocity(istringstream &stream, vector<Settings> &settings){
	vector<double> values;
	string value_string = "";
	int flag = 0;
	long unsigned int n_settings = settings.size();
	long unsigned int n_values = 0;

	// Read values for the target velocity
	getline(stream, value_string, DELIMITER);
	flag = readDoubles(values, value_string);

	if(flag == 0){
		n_settings = settings.size();
		for(unsigned int i = 0; i < n_settings; ++i){
			settings[i].velocity.push_back(values[0]);
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
		}
		for(unsigned int i = 0; i < values.size(); ++i){
			for(unsigned int j = 0; j < n_settings; ++j){
				settings[i*n_settings + j].velocity.push_back(values[i]);
			}
		}
	}
	if(flag == 2){
		n_values = values.size();

		if(n_settings == 1){
			for(unsigned int i = 0; i < n_values - 1; ++i){
				settings.push_back(settings[0]);
			}
		}

		for(unsigned int i = 0; i < n_values; ++i){
			settings[i].velocity.push_back(values[i]);
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
        //cout << "> Reading input file '" << filename.str() << "'" << endl;

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

	stringstream filename;
	filename << MU_DIR << mass_attenuation_filename;
	ifstream ifile(filename.str());

        if(!ifile.is_open()){
		cout << "Error: " << __FILE__ << ":" << __LINE__ << ": "; 
		cout << " readNIST(): File '" << filename.str() << "' not found." << endl;
		abort();
	}
        //cout << "> Reading input file '" << mass_attenuation_filename << "'" << endl;

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

void InputReader::read1ColumnFile(vector<double> &data, string filename){

	string line, value;
	ifstream ifile;

	ifile.open(filename);	

        if(!ifile.is_open()){
		cout << "Error: " << __FILE__ << ":" << __LINE__ << ": "; 
		cout << " read2ColumnFile(): File '" << filename << "' not found." << endl;
		abort();
	}
        cout << "> Reading input file '" << filename << "'" << endl;

	while(getline(ifile, line)){
		if(line.substr(0,1) == COMMENT)
			continue;

		istringstream stream(line);

		getline(stream, value, DELIMITER);
		data.push_back(atof(value.c_str()));
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
        //cout << "> Reading input file '" << filename << "'" << endl;

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

void InputReader::read3ColumnFile(vector< vector<double> > &data, string filename){

	string line, value;
	ifstream ifile;

	ifile.open(filename);	

        if(!ifile.is_open()){
		cout << "Error: " << __FILE__ << ":" << __LINE__ << ": "; 
		cout << " read3ColumnFile(): File '" << filename << "' not found." << endl;
		abort();
	}
        //cout << "> Reading input file '" << filename << "'" << endl;

	data.push_back(vector<double>());
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
		getline(stream, value, DELIMITER); 
		data[2].push_back(atof(value.c_str()));
	}
}
