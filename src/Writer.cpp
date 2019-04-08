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


#include "Config.h"
#include "Writer.h"

#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>

using std::stringstream;
using std::ofstream;
using std::cout;
using std::endl;
using std::scientific;
using std::defaultfloat;
using std::setprecision;

void Writer::write1DHistogram(const vector<double> &histogram, const string name, const string column_name){

	stringstream filename;
	filename << settings.outputfile << "_" << name << TXT_SUFFIX;

	ofstream ofile(filename.str());

        if(!ofile.is_open()){
                cout << "Error: " << __FILE__ << ":" << __LINE__ << ": "; 
		cout << " write1DHistogram(): File '" << filename.str() << "' could not be opened." << endl;
		abort();
	}
        cout << "> Writing output file '" << filename.str() << "'" << endl;

	ofile.precision(8);
	ofile << COMMENT << " " << column_name << endl;
	for(unsigned int i = 0; i < histogram.size(); ++i){
		ofile << scientific << histogram[i] << endl;
	}
}

void Writer::write2DHistogram(const vector<vector<double> > &histogram, const string name, const string line_name, const string column_name){

	stringstream filename;
	filename << settings.outputfile << "_" << name << TXT_SUFFIX;

	ofstream ofile(filename.str());

        if(!ofile.is_open()){
                cout << "Error: " << __FILE__ << ":" << __LINE__ << ": "; 
		cout << " write2DHistogram(): File '" << filename.str() << "' could not be opened." << endl;
		abort();
	}
        cout << "> Writing output file '" << filename.str() << "'" << endl;

	ofile.precision(8);
	ofile << COMMENT << " Lines   : " << line_name << endl;
	ofile << COMMENT << " Columns : " << column_name << endl;

	long unsigned int ncolumns = histogram.size();
	long unsigned int nlines = histogram[0].size();

	for(unsigned int i = 0; i < ncolumns; ++i){
		for(unsigned int j = 0; j < nlines; ++j){
			ofile << scientific << histogram[i][j] << "\t";
		}
		ofile << endl;
	}
}

void Writer::write1DCalibration(const vector<double> &bins, const string name, const string histogram_name){
	unsigned int nbins = (unsigned int) bins.size();
	double a = (bins[nbins - 1] - bins[0])/((double) nbins - 1.);
	double b = bins[0];

	stringstream filename;
	filename << settings.outputfile << "_" << name << CAL_SUFFIX;

	ofstream ofile(filename.str(), std::ios_base::app);

        if(!ofile.is_open()){
                cout << "Error: " << __FILE__ << ":" << __LINE__ << ": "; 
		cout << " write2DHistogram(): File '" << filename.str() << "' could not be opened." << endl;
		abort();
	}
        cout << "> Writing output file '" << filename.str() << "'" << endl;

	ofile.precision(8);
	ofile << scientific << histogram_name << TXT_SUFFIX << ":\t" << b << "\t" << a << endl;
	ofile.close();
}
