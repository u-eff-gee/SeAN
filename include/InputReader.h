#ifndef INPUTREADER_H
#define INPUTREADER_H 1

#include "Settings.h"

#include <string>

using std::string;

class InputReader{
public:
	InputReader(){};
	~InputReader(){};

	// Read the SeAN input file
	void readFile(Settings &settings);
	// Read a nuclear mass from the AME table
	double readAME(string isotope);
	// Read an energy-dependent mass-attenuation coefficient from a file with the NIST format
	void readNIST(vector< vector<double> > &mass_attenuation_file, const string mass_attenuation_filename);
	// Read any file with two columns
	void read2ColumnFile(vector< vector<double> > &data, string filename);


private:
	// Settings for the SeAN input file
	unsigned int N_EXPERIMENT_SETTINGS = 4; 
	unsigned int N_TARGET_SETTINGS = 11;
	char DELIMITER = ',';

	// Settings for the directories
	string V_DIR = "velocity_distribution/";
	string MU_DIR = "mass_attenuation/";
	string BEAM_DIR = "beam/";
	string MASS_DIR = "atomic_mass/";
	string AME_FILE_NAME = "mass16.txt";

	// Coordinates in AME file
const 	unsigned int AME_HEADER_LENGTH = 38;
const 	unsigned int AME_MASS_NUMBER = 16;
const 	unsigned int AME_ISOTOPE = 20;
const 	unsigned int AME_MASS_START = 96;
const 	unsigned int AME_MASS_LENGTH = 13;
	
	// Coordinates in NIST file
const 	unsigned int NIST_XRAY =  0;
const 	unsigned int NIST_XRAY_LENGTH = 3;
const 	unsigned int NIST_ENERGY = 3;
const 	unsigned int NIST_ENERGY_LENGTH = 11;
const 	unsigned int NIST_MU = 16;
const 	unsigned int NIST_MU_LENGTH = 9;

};

#endif
