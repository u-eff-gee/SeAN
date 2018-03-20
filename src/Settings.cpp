#include "Settings.h"
#include "Config.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

using std::cout;
using std::endl;
using std::scientific;
using std::fixed;
using std::defaultfloat;
using std::setprecision;
using std::stringstream;
using std::ofstream;

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
	if(inputfile != ""){
		cout << "FILE:\t" << inputfile << endl;
	} else{
		cout << "FILE:\tnot set" << endl;
	}
	cout << "COMMAND LINE OPTIONS: " << endl;
	cout << "\tEXACT:\t" << exact << endl;
	cout << "\tPLOT :\t" << plot << endl;
	if(write){
		cout << "\tWRITE:\t" << "true" << endl;
	}
	else{
		cout << "\tWRITE:\t" << "false" << endl;
	}
	cout << "\tVERBOSITY:\t" << verbosity << endl;
	if(output){
		if(outputfile != ""){
			cout << "\tOUTPUT:\t" << outputfile << endl;
		}
		else{
			cout << "\tOUTPUT:\ttrue" << endl;
		}
	}
	else
		cout << "\tOUTPUT:\t" << "false" << endl;
}

void Settings::printExperiment(){
	long int default_precision = cout.precision();

	cout << HORIZONTAL_LINE << endl;
	cout << ">>> EXPERIMENT" << endl;
	// Since eV is the standard energy unit of SeAN, but the absolute energy of resonances is in the order of several MeV, increase the precision for the output of double numbers and reset it later
	cout << setprecision(7) << "EMIN\t\t:\t" << scientific << emin << " eV" << endl;
	cout << setprecision(7) << "EMAX\t\t:\t" << scientific << emax << " eV" << endl;
	cout << "PRIMARY BEAM\t:\t";
	
	// Reset precision
	cout << defaultfloat << setprecision((int) default_precision);

	if(incidentBeamParams.size()){
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
	} else{
		cout << "not set" << endl;
	}

	cout << "NBINS_ENERGY\t:\t" << nbins_e << endl;
	cout << "NBINS_Z\t\t:\t" << nbins_z << endl;

}

void Settings::printTarget(unsigned int i){
	long int default_precision = cout.precision();

	cout << HORIZONTAL_LINE << endl;
	cout << ">>> TARGET #" << (i + 1) << " : " << targetNames[i]  << endl;
	cout << "RESONANCES:\tENERGY\tGAMMA0\tGAMMA\tJ0\tJ" << endl;

	if(energy.size() == 0 || gamma0.size() == 0 || gamma.size() == 0 || jj.size() == 0){
		cout << "\tno resonances given or incomplete input" << endl;
	} else{

		long unsigned int nresonances = energy[i].size();

		if(energy[i].size() == nresonances && gamma0[i].size() == nresonances && gamma[i].size() && jj[i].size() == nresonances){

			for(long unsigned int j = 0; j < nresonances; ++j){
				cout << "\t\t" << scientific << setprecision(7) << energy[i][j] << defaultfloat << "\t" << gamma0[i][j] << "\t" << gamma[i][j] << "\t" << ji[i] << "\t" << jj[i][j] << endl;
			}
			
			// Reset precision
			cout << defaultfloat << setprecision((int) default_precision);
		} else{
			cout << "Error: " << __FILE__ << ":" << __LINE__ << ": "; 
			cout << " printTarget(): Number of parameters ENERGY, GAMMA0, GAMMA, J0 and J does not match." << endl;
			abort();
		}
	}

	cout << "DOPPLER_BROADENING:\t";

	if(dopplerBroadening.size() == 0){
		cout << "not set" << endl;
	} else{
		switch(dopplerBroadening[i]){
			case dopplerModel::zero:
				cout << "zero" << endl;
				break;
			case dopplerModel::arb_vdist:
				cout << "arb_vdist, " << vDistFile[i] << endl;
				break;
			case dopplerModel::arb_cs:
				cout << "arb_cs, " << dopplerFile[i] << endl;
				break;
			case dopplerModel::mb:
				cout << "Maxwell-Boltzmann, T_eff = " << dopplerParams[i][0] << " K" << endl;
				break;
			case dopplerModel::mba:
				cout << "Maxwell-Boltzmann (using approximation), T_eff = " << dopplerParams[i][0] << " K " << endl;
				break;
			case dopplerModel::mbd:
				cout << "Maxwell-Boltzmann (using Debye approximation), T = " << dopplerParams[i][0] << " K, T_D = " << dopplerParams[i][1] << " K " << endl;
				break;
			case dopplerModel::mbad:
				cout << "Maxwell-Boltzmann (using Debye and integral approximation), T = " << dopplerParams[i][0] << " K, T_D = " << dopplerParams[i][1] << " K " << endl;
				break;
			default: break;
		}
	}
	
	cout << "ATOMIC MASS:\t\t";
	if(mass.size() == 0){
		cout << "not set" << endl;
	} else{
		cout << mass[i] << " u" << endl;
	}
	
	cout << "MASS ATTENUATION:\t";

	if(mAttParams.size() == 0 && mAttFile.size() == 0){
		cout << "not set" << endl;
	} else{
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
	}
	
	cout << "TARGET THICKNESS:\t";
	if(thickness.size() == 0){
		cout << "not set" << endl;
	} else{
	 	cout << thickness[i] << " atoms/fm^2" << endl;
	}

	cout << "VELOCITY:\t\t";
	if(velocity.size() == 0){
		cout << "not set" << endl;
	} else{
		cout << velocity[i] << " m/s" << endl;
	}
}

void Settings::write_output(unsigned int n_setting) const {

	stringstream filename;
	filename << "output/" << outputfile;

	ofstream ofile;
	ofile.open(filename.str(), std::ios_base::out | std::ios_base::app);

        if(!ofile.is_open()){
                cout << "Error: " << __FILE__ << ":" << __LINE__ << ": "; 
		cout << " write_output(): File '" << outputfile << "' could not be opened." << endl;
		abort();
	}

	writeOptions(n_setting);
	writeExperiment();
	unsigned int ntargets = (unsigned int) targetNames.size();
	for(unsigned int i = 0; i < ntargets; ++i)
		writeTarget(i);

	ofile << HORIZONTAL_LINE << endl;

	ofile.close();
}

void Settings::writeOptions(unsigned int n_setting) const{

	stringstream filename;
	filename << "output/" << outputfile;

	ofstream ofile;
	ofile.open(filename.str(), std::ios_base::out | std::ios_base::app);

        if(!ofile.is_open()){
                cout << "Error: " << __FILE__ << ":" << __LINE__ << ": "; 
		cout << " write_input(): File '" << outputfile << "' could not be opened." << endl;
		abort();
	}

	ofile << HORIZONTAL_LINE << endl;
	ofile << ">>> SeAN INPUT #" << n_setting << endl;
	if(inputfile != ""){
		ofile << "FILE:\t" << inputfile << endl;
	} else{
		ofile << "FILE:\tnot set" << endl;
	}
	ofile << "COMMAND LINE OPTIONS: " << endl;
	ofile << "\tEXACT:\t" << exact << endl;
	ofile << "\tPLOT :\t" << plot << endl;
	if(write){
		ofile << "\tWRITE:\t" << "true" << endl;
	}
	else{
		ofile << "\tWRITE:\t" << "false" << endl;
	}
	ofile << "\tVERBOSITY:\t" << verbosity << endl;
	if(output){
		if(outputfile != ""){
			ofile << "\tOUTPUT:\t" << outputfile << endl;
		}
		else{
			ofile << "\tOUTPUT:\ttrue" << endl;
		}
	}
	else
		ofile << "\tOUTPUT:\t" << "false" << endl;

	ofile.close();
}

void Settings::writeExperiment() const{

	stringstream filename; filename << "output/" << outputfile;

	ofstream ofile;
	ofile.open(filename.str(), std::ios_base::out | std::ios_base::app);

        if(!ofile.is_open()){
                cout << "Error: " << __FILE__ << ":" << __LINE__ << ": "; 
		cout << " write_input(): File '" << outputfile << "' could not be opened." << endl;
		abort();
	}

	long int default_precision = ofile.precision();

	ofile << HORIZONTAL_LINE << endl;
	ofile << ">>> EXPERIMENT" << endl;
	// Since eV is the standard energy unit of SeAN, but the absolute energy of resonances is in the order of several MeV, increase the precision for the output of double numbers and reset it later
	ofile << setprecision(7) << "EMIN\t\t:\t" << scientific << emin << " eV" << endl;
	ofile << setprecision(7) << "EMAX\t\t:\t" << scientific << emax << " eV" << endl;
	ofile << "PRIMARY BEAM\t:\t";
	
	// Reset precision
	ofile << defaultfloat << setprecision((int) default_precision);

	if(incidentBeamParams.size()){
		switch(incidentBeam){
			case incidentBeamModel::constant:
				ofile << "constant, " << incidentBeamParams[0] << endl;
				break;
			case incidentBeamModel::gauss:
				ofile << "gauss, MU = " << incidentBeamParams[0] << " eV, SIGMA = " << incidentBeamParams[1] << " eV, NORM = " << incidentBeamParams[2] << endl;
				break;
			case incidentBeamModel::arb:
				ofile << "arb, " << incidentBeamFile << endl;
				break;
			default: break;
		}
	} else{
		ofile << "not set" << endl;
	}

	ofile << "NBINS_ENERGY\t:\t" << nbins_e << endl;
	ofile << "NBINS_Z\t\t:\t" << nbins_z << endl;
}

void Settings::writeTarget(unsigned int i) const{

	stringstream filename; filename << "output/" << outputfile;

	ofstream ofile;
	ofile.open(filename.str(), std::ios_base::out | std::ios_base::app);

        if(!ofile.is_open()){
                cout << "Error: " << __FILE__ << ":" << __LINE__ << ": "; 
		cout << " write_input(): File '" << outputfile << "' could not be opened." << endl;
		abort();
	}

	long int default_precision = ofile.precision();

	ofile << HORIZONTAL_LINE << endl;
	ofile << ">>> TARGET #" << (i + 1) << " : " << targetNames[i]  << endl;
	ofile << "RESONANCES:\tENERGY\tGAMMA0\tGAMMA\tJ0\tJ" << endl;

	if(energy.size() == 0 || gamma0.size() == 0 || gamma.size() == 0 || jj.size() == 0){
		ofile << "\tno resonances given or incomplete input" << endl;
	} else{

		long unsigned int nresonances = energy[i].size();

		if(energy[i].size() == nresonances && gamma0[i].size() == nresonances && gamma[i].size() && jj[i].size() == nresonances){

			for(long unsigned int j = 0; j < nresonances; ++j){
				ofile << "\t\t" << scientific << setprecision(7) << energy[i][j] << defaultfloat << "\t" << gamma0[i][j] << "\t" << gamma[i][j] << "\t" << ji[i] << "\t" << jj[i][j] << endl;
			}
			
			// Reset precision
			ofile << defaultfloat << setprecision((int) default_precision);
		} else{
			ofile << "Error: " << __FILE__ << ":" << __LINE__ << ": "; 
			ofile << " printTarget(): Number of parameters ENERGY, GAMMA0, GAMMA, J0 and J does not match." << endl;
			abort();
		}
	}

	ofile << "DOPPLER BROADENING:\t";

	if(dopplerBroadening.size() == 0){
		ofile << "not set" << endl;
	} else{
		switch(dopplerBroadening[i]){
			case dopplerModel::zero:
				ofile << "zero" << endl;
				break;
			case dopplerModel::arb_vdist:
				ofile << "arb_vdist, " << dopplerFile[i] << endl;
				break;
			case dopplerModel::arb_cs:
				ofile << "arb_cs, " << dopplerFile[i] << endl;
				break;
			case dopplerModel::mb:
				ofile << "Maxwell-Boltzmann, T_eff = " << dopplerParams[i][0] << " K" << endl;
				break;
			case dopplerModel::mba:
				ofile << "Maxwell-Boltzmann (using approximation), T_eff = " << dopplerParams[i][0] << " K " << endl;
				break;
			case dopplerModel::mbd:
				ofile << "Maxwell-Boltzmann (using Debye approximation), T = " << dopplerParams[i][0] << " K, T_D = " << dopplerParams[i][1] << " K " << endl;
				break;
			case dopplerModel::mbad:
				ofile << "Maxwell-Boltzmann (using Debye and integral approximation), T = " << dopplerParams[i][0] << " K, T_D = " << dopplerParams[i][1] << " K " << endl;
				break;
			default: break;
		}
	}
	
	ofile << "ATOMIC MASS:\t\t";
	if(mass.size() == 0){
		ofile << "not set" << endl;
	} else{
		ofile << mass[i] << " u" << endl;
	}
	
	ofile << "MASS ATTENUATION:\t";

	if(mAttParams.size() == 0 && mAttFile.size() == 0){
		ofile << "not set" << endl;
	} else{
		switch(mAtt[i]){
			case mAttModel::constant:
				ofile << "constant, " << mAttParams[i][0] << endl;
				break;
			case mAttModel::nist:
				ofile << "nist, " << mAttFile[i] << endl;
				break;
			case mAttModel::arb:
				ofile << "arb" << mAttFile[i] << endl;
				break;
			default: break;
		}
	}
	
	ofile << "TARGET THICKNESS:\t";
	if(thickness.size() == 0){
		ofile << "not set" << endl;
	} else{
	 	ofile << thickness[i] << " atoms/fm^2" << endl;
	}

	ofile << "VELOCITY:\t\t";
	if(velocity.size() == 0){
		ofile << "not set" << endl;
	} else{
		ofile << velocity[i] << " m/s" << endl;
	}
}
