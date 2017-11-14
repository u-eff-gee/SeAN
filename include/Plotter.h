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

	void plot1DHistogram(vector<double> &bins, vector<double> &histogram, string name);
};

#endif
