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
};

#endif
