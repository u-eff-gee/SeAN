#include "Settings.h"

#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::scientific;
using std::setprecision;

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
	cout << "###################################################" << endl;
	cout << ">>> EXPERIMENT" << endl;
	cout << setprecision(7) << "EMIN\t\t:\t" << scientific << emin << " eV" << endl;
	cout << setprecision(7) << "EMAX\t\t:\t" << scientific << emax << " eV" << endl;
	cout << "NBINS ENERGY\t:\t" << nbins_e << endl;
	cout << "NBINS Z\t\t:\t" << nbins_z << endl;
	cout << "PRIMARY BEAM\t:\t";
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
