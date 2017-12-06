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

	void plot1DHistogram(const vector<double> &bins, const vector<double> &histogram, const string name, const string xaxis_label, const string yaxis_label);

	void plotMultiple1DHistograms(const vector<double> &bins, const vector< vector<double> > &histograms, const string name, const string xaxis_label, const string yaxis_label);

	void plotMultiple1DHistogramsAndSum(const vector<double> &bins, const vector< vector<double> > &histograms, const string name, const string xaxis_label, const string yaxis_label);

	void plot2DHistogram(const vector<double> &bins1, const vector<double> &bins2, const vector< vector<double> > &histogram, const string name, const string xaxis_label, const string yaxis_label, const string zaxis_label);
};

#endif
