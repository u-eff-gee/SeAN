#ifndef WRITER_H
#define WRITER_H 1

#include <string>
#include <vector>

using std::string;
using std::vector;

class Writer{

private:
	string TXT_OUTPUT_DIR = "output/";
	const string COMMENT = "#";

public:
	Writer(){};
	~Writer(){};

	void write1DHistogram(const vector<double> &histogram, const string name, const string column_name);

	void write2DHistogram(const vector<vector<double> > &histogram, const string name, const string column1_name, const string column2_name);
};

#endif
