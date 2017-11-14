#ifndef WRITER_H
#define WRITER_H 1

#include <string>
#include <vector>

using std::string;
using std::vector;

class Writer{

private:
	string TXT_OUTPUT_DIR = "output/";

public:
	Writer(){};
	~Writer(){};

	void write1DHistogram(vector<double> &histogram, string name, string column_name);
};

#endif
