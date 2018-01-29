#ifndef INPUTREADER_H
#define INPUTREADER_H 1

#include "Settings.h"

#include <string>

using std::string;
using std::istringstream;

class InputReader{
public:
	InputReader(){};
	~InputReader(){};

	// Read the SeAN input file
	void readFile(vector<Settings> &settings);
	// Reader functions for parameters
		// General reader function for double values
	int readDoubles(vector<double> &values, const string &value_string);
		// General reader function for int values
	int readInts(vector<int> &values, const string &value_string);

	void readEminEmax(istringstream &stream, vector<Settings> &settings);
	void readIncidentBeam(istringstream &stream, vector<Settings> &settings);
	void readNBins_E(istringstream &stream, vector<Settings> &settings);
	void readNBins_Z(istringstream &stream, vector<Settings> &settings);
	void readEnergy(istringstream &stream, vector<Settings> &settings, const unsigned int ntarget);

	void readGamma0(istringstream &stream, vector<Settings> &settings, const unsigned int ntarget);

	void readGamma(istringstream &stream, vector<Settings> &settings, const unsigned int ntarget);

	void readJ0(istringstream &stream, vector<Settings> &settings, const unsigned int ntarget);

	void readJ(istringstream &stream, vector<Settings> &settings, const unsigned int ntarget);

	// Read a nuclear mass from the AME table
	double readAME(string isotope);
	// Read an energy-dependent mass-attenuation coefficient from a file with the NIST format
	void readNIST(vector< vector<double> > &mass_attenuation_file, const string mass_attenuation_filename);
	// Read any file with two columns
	void read2ColumnFile(vector< vector<double> > &data, string filename);


private:
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
