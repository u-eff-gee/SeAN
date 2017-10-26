#ifndef INPUTREADER_H
#define INPUTREADER_H 1

#include "Settings.h"

#include <string>

using std::string;

class InputReader{
public:
	InputReader(){};
	~InputReader(){};

	void readFile(const char* filename, Settings &settings);

private:
	unsigned int N_EXPERIMENT_SETTINGS = 4; 
	unsigned int N_TARGET_SETTINGS = 11;
	char DELIMITER = ',';
};

#endif
