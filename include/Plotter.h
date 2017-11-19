#ifndef PLOTTER_H
#define PLOTTER_H 1

#include <string>
#include <vector>

using std::string;
using std::vector;

class Plotter{
private:
	string PLOT_OUTPUT_DIR = "output/";

public:
	Plotter(){};
	~Plotter(){};

	void plot1DHistogram(const vector<double> &bins, const vector<double> &histogram, const string name);

	void plotMultiple1DHistograms(const vector<double> &bins, const vector< vector<double> > &histograms, const string name);

	void plotMultiple1DHistogramsAndSum(const vector<double> &bins, const vector< vector<double> > &histograms, string name);
};

#endif
