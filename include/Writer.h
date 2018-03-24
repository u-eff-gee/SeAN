#ifndef WRITER_H
#define WRITER_H 1

#include <string>
#include <vector>

#include "Settings.h"

using std::string;
using std::vector;

class Writer{

public:
	Writer(){};
	~Writer(){};

	// Methods to write histograms to a file
	void write1DHistogram(const vector<double> &histogram, const string name, const string column_name);

	void write2DHistogram(const vector<vector<double> > &histogram, const string name, const string column1_name, const string column2_name);

	// Methods to write calibration parameters for the histograms that allow the user to convert a bin number into a physical quantity. For example, in a cross section histogram, if bin #0 corresponds to 10 eV and bin #1 corresponds to 30 eV, the calibration parameters would be a*#bin + b with a = 20. and b = 10.
	void write1DCalibration(const vector<double> &bins, const string name);
	// Find out how to append to a file
};

#endif
