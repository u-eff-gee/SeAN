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
	filename << OUTPUT_DIR << name << ".txt";

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
	filename << OUTPUT_DIR << name << ".txt";

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

	for(unsigned int i = 0; i < histogram.size(); ++i){
		for(unsigned int j = 0; j < histogram.size(); ++j){
			ofile << scientific << histogram[i][j] << ", ";
		}
		ofile << endl;
	}
}
