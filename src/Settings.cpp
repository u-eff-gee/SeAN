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

	cout << HORIZONTAL_LINE << endl;
}

void Settings::printOptions(){
	cout << HORIZONTAL_LINE << endl;
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

	cout << HORIZONTAL_LINE << endl;
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
			cout << "arb, " << incidentBeamFile << endl;
			break;
		default: break;
	}
}

void Settings::printTarget(unsigned int i){
	long int default_precision = cout.precision();

	cout << HORIZONTAL_LINE << endl;
	cout << ">>> TARGET #" << (i + 1) << " : " << targetNames[i]  << endl;
	cout << "RESONANCES:\tENERGY\tGAMMA0\tGAMMA\tJ0\tJ" << endl;

	long unsigned int nresonances = energy[i].size();
	for(long unsigned int j = 0; j < nresonances; ++j){
		cout << "\t\t" << scientific << setprecision(7) << energy[i][j] << defaultfloat << "\t" << gamma0[i][j] << "\t" << gamma[i][j] << "\t" << ji[i] << "\t" << jj[i][j] << endl;
	}
	
	// Reset precision
	cout << defaultfloat << setprecision((int) default_precision);

	cout << "VELOCITY DISTRIBUTION:\t";
	switch(vDist[i]){
		case vDistModel::zero:
			cout << "zero" << endl;
			break;
		case vDistModel::arb:
			cout << "arb, " << vDistFile[i] << endl;
			break;
		case vDistModel::mb:
			cout << "Maxwell-Boltzmann, T_eff = " << vDistParams[i][0] << " K" << endl;
			break;
		case vDistModel::mba:
			cout << "Maxwell-Boltzmann (using approximation), T_eff = " << vDistParams[i][0] << " K " << endl;
			break;
		case vDistModel::mbd:
			cout << "Maxwell-Boltzmann (using Debye approximation), T = " << vDistParams[i][0] << " K, T_D = " << vDistParams[i][1] << " K " << endl;
			break;
		case vDistModel::mbad:
			cout << "Maxwell-Boltzmann (using Debye and integral approximation), T = " << vDistParams[i][0] << " K, T_D = " << vDistParams[i][1] << " K " << endl;
			break;
		default: break;
	}

	cout << "ATOMIC MASS:\t\t" << mass[i] << " u" << endl;

	cout << "MASS ATTENUATION:\t";
	switch(mAtt[i]){
		case mAttModel::constant:
			cout << "constant, " << mAttParams[i][0] << endl;
			break;
		case mAttModel::nist:
			cout << "nist, " << mAttFile[i] << endl;
			break;
		case mAttModel::arb:
			cout << "arb" << mAttFile[i] << endl;
			break;
		default: break;
	}

	cout << "TARGET THICKNESS:\t" << thickness[i] << " atoms/fm^2" << endl;

	cout << "VELOCITY:\t\t" << velocity[i] << " m/s" << endl;
}
