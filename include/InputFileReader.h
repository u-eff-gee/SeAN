#ifndef INPUTFILEREADER_H 
#define INPUTFILEREADER_H 1

#include <vector>

#include "Target.h"

using std::vector;

class InputFileReader{

private:

	vector<Target*> targets;
	vector<double> beamParams;

	string ifname;
	string beam_ID;

	double emin;
	double emax;

public:
	InputFileReader(){};
	~InputFileReader(){};

	void readInputFile(const char* filename);
	void print();
};

#endif
