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

	// Coordinates in AME file
	unsigned int AME_HEADER_LENGTH = 38;
	unsigned int AME_MASS_NUMBER = 16;
	unsigned int AME_ISOTOPE = 20;
	unsigned int AME_MASS_START = 96;
	unsigned int AME_MASS_LENGTH = 13;
};

#endif
