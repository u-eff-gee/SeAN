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
};

#endif
