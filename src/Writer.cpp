#include "Config.h"
#include "Writer.h"

#include <sstream>
#include <fstream>
#include <iostream>

using std::stringstream;
using std::ofstream;
using std::cout;
using std::endl;
using std::scientific;

void Writer::write1DHistogram(vector<double> &histogram, string name, string column_name){

	stringstream filename;
	filename << TXT_OUTPUT_DIR << name << ".txt";

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
