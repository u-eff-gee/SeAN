#include "Settings.h"

#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::scientific;
using std::fixed;
using std::defaultfloat;
using std::setprecision;

void Settings::print(){

	printOptions();
	printExperiment();
	unsigned int ntargets = (unsigned int) targetNames.size();
	for(unsigned int i = 0; i < ntargets; ++i)
		printTarget(i);
}

void Settings::printOptions(){
	cout << "###################################################" << endl;
	cout << ">>> SeAN INPUT" << endl;
	cout << "FILE:\t" << inputfile << endl;
	cout << "COMMAND LINE OPTIONS: " << endl;
	cout << "\tEXACT:\t" << exact << endl;
	cout << "\tPLOT :\t" << plot << endl;
	if(write || sudowrite){
		cout << "\tWRITE:\t" << "1" << endl;
	}
	else{
		cout << "\tWRITE:\t" << "0" << endl;
	}
}

void Settings::printExperiment(){
	long int default_precision = cout.precision();

	cout << "###################################################" << endl;
	cout << ">>> EXPERIMENT" << endl;
	// Since eV is the standard energy unit of SeAN, but the absolute energy of resonances is in the order of several MeV, increase the precision for the output of double numbers and reset it later
	cout << setprecision(7) << "EMIN\t\t:\t" << scientific << emin << " eV" << endl;
	cout << setprecision(7) << "EMAX\t\t:\t" << scientific << emax << " eV" << endl;
	cout << "NBINS_ENERGY\t:\t" << nbins_e << endl;
	cout << "NBINS_Z\t\t:\t" << nbins_z << endl;
	cout << "PRIMARY BEAM\t:\t";

	// Reset precision
	cout << defaultfloat << setprecision((int) default_precision);

	switch(incidentBeam){
		case incidentBeamModel::constant:
			cout << "constant, " << incidentBeamParams[0] << endl;
			break;
		case incidentBeamModel::gauss:
			cout << "gauss, MU = " << incidentBeamParams[0] << " eV, SIGMA = " << incidentBeamParams[1] << " eV, NORM = " << incidentBeamParams[2] << endl;
			break;
		case incidentBeamModel::arb:
			break;
		default: break;
	}
}

void Settings::printTarget(unsigned int i){
	long int default_precision = cout.precision();

	cout << "###################################################" << endl;
	cout << ">>> TARGET #" << (i + 1) << " : " << targetNames[i]  << endl;
	cout << "RESONANCES:\tENERGY\tGAMMA0\tGAMMA\tJ0\tJ" << endl;

	long unsigned int nresonances = energy[i].size();
	for(long unsigned int j = 0; j < nresonances; ++j){
		cout << "\t\t" << scientific << setprecision(7) << energy[i][j] << defaultfloat << "\t" << gamma0[i][j] << "\t" << gamma[i][j] << "\t" << ji[i] << "\t" << jj[i][j] << endl;
	}
	
	// Reset precision
	cout << defaultfloat << setprecision((int) default_precision);

	cout << "MASS ATTENUATION\t:\t";
	switch(mAtt){
		case mAttModel::constant:
			cout << "constant, " << mAttParams[0] << endl;
			break;
		case mAttModel::nist:
			cout << "nist, " << mAttFile << endl;
			break;
		case mAttModel::arb:
			break;
		default: break;
	}

}
