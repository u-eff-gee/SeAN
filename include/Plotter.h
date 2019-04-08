/*    
    This file is part of SeAN.

    SeAN is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SeAN is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SeAN.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef PLOTTER_H
#define PLOTTER_H 1

#include <string>
#include <vector>

#include "Settings.h"

using std::string;
using std::vector;

class Plotter{

private:
	Settings settings;

public:
	Plotter(Settings &s){
		settings = s;
	};
	~Plotter(){};

	void plot1DHistogram(const vector<double> &bins, const vector<double> &histogram, const string name, const string xaxis_label, const string yaxis_label);

	void plotMultiple1DHistograms(const vector<double> &bins, const vector< vector<double> > &histograms, const string name, const string xaxis_label, const string yaxis_label);

	void plotMultiple1DHistogramsAndSum(const vector<double> &bins, const vector< vector<double> > &histograms, const string name, const string xaxis_label, const string yaxis_label);

	void plot2DHistogram(const vector<double> &bins1, const vector<double> &bins2, const vector< vector<double> > &histogram, const string name, const string xaxis_label, const string yaxis_label, const string zaxis_label);
};

#endif
