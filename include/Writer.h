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

#pragma once 

#include <string>
#include <vector>

#include "Settings.h"

using std::string;
using std::vector;

class Writer{
private:
	Settings settings;

public:
	Writer(Settings &s){
		settings = s;
	};
	~Writer(){};

	// Methods to write histograms to a file
	void write1DHistogram(const vector<double> &histogram, const string name, const string column_name);

	void write2DHistogram(const vector<vector<double> > &histogram, const string name, const string column1_name, const string column2_name);

	// Methods to write calibration parameters for the histograms that allow the user to convert a bin number into a physical quantity. For example, in a cross section histogram, if bin #0 corresponds to 10 eV and bin #1 corresponds to 30 eV, the calibration parameters would be a*#bin + b with a = 20. and b = 10.
	void write1DCalibration(const vector<double> &bins, const string name, const string histogram_name);
	// Find out how to append to a file
};