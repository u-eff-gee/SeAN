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

string Settings::bool_string(const bool b) const {
		if(b) return "true"; else return "false";
}

void Settings::print(){

	printOptions();
	printExperiment();
	unsigned int ntargets = (unsigned int) targetNames.size();
	for(unsigned int i = 0; i < ntargets; ++i)
		printTarget(i);
}

string Settings::option_string() const {

	stringstream opss;
	
	opss << HORIZONTAL_LINE << "\n";
	opss << ">>> SeAN INPUT\n";
	if(inputfile != ""){
		opss << "FILE:\t" << inputfile << "\n" ;
	} else{
		opss << "FILE:\tnot set\n" ;
	}
	opss << "COMMAND LINE OPTIONS: \n" ;
	opss << "\tDIRECT      :\t" << direct << "\n";
	opss << "\tEXACT       :\t" << exact << "\n";
	opss << "\tPLOT        :\t" << plot << "\n";
	opss << "\tMULTI       :\t" << bool_string(multi) << "\n";
	opss << "\tWRITE       :\t" << bool_string(write);
	if(write_all) opss << " ( ALL )" << "\n"; else opss << "\n";
	opss << "\tRECOIL      :\t" << bool_string(recoil) << "\n" ;
	opss << "\tUNCERTAINTY :\t" << bool_string(uncertainty) << "\n" ;
	opss << "\tVERBOSITY   :\t" << verbosity << "\n";
	opss << "\tOUTPUT      :\t" << bool_string(output);
	if(output) opss << " ( " << outputfile << " )";
	opss << "\n"; 

	return opss.str();
};

string Settings::experiment_string() const {
	long int default_precision = cout.precision();

	stringstream exss;

	exss << HORIZONTAL_LINE << "\n";
	exss << ">>> EXPERIMENT\n" ;
	// Since eV is the standard energy unit of SeAN, but the absolute energy of resonances is in the order of several MeV, increase the precision for the output of double numbers and reset it later
	exss << setprecision(7) << "EMIN\t\t:\t" << scientific << emin << " eV\n" ;
	exss << setprecision(7) << "EMAX\t\t:\t" << scientific << emax << " eV\n" ;
	exss << "PRIMARY BEAM\t:\t";
	
	// Reset precision
	exss << defaultfloat << setprecision((int) default_precision);

	if(incidentBeamParams.size()){
		switch(incidentBeam){
			case incidentBeamModel::constant:
				exss << "constant, " << incidentBeamParams[0] ;
				break;
			case incidentBeamModel::gauss:
				exss << "gauss, MU = " << incidentBeamParams[0] << " eV, SIGMA = " << incidentBeamParams[1] << " eV, NORM = " << incidentBeamParams[2] ;
				break;
			case incidentBeamModel::arb:
				exss << "arb, " << incidentBeamFile ;
				break;
			default: break;
		}
	} else{
		exss << "not set" ;
	}

	exss << "\nNBINS_ENERGY\t:\t" << nbins_e ;
	exss << "\nNBINS_Z\t\t:\t" << nbins_z <<"\n";

	return exss.str();
}

string Settings::target_string(unsigned int i) const {

	long int default_precision = cout.precision();

	stringstream tass;

	tass << HORIZONTAL_LINE << "\n";
	tass << ">>> TARGET #" << (i + 1) << " : " << targetNames[i] << "\n" ;
	tass << "RESONANCES:\tENERGY\tGAMMA0\tGAMMA\tJ0\tJ\n" ;

	if(energy.size() == 0 || gamma0.size() == 0 || gamma.size() == 0 || jj.size() == 0){
		tass << "\tno resonances given or incomplete input\n" ;
	} else{

		long unsigned int nresonances = energy[i].size();

		if(energy[i].size() == nresonances && gamma0[i].size() == nresonances && gamma[i].size() && jj[i].size() == nresonances){

			for(long unsigned int j = 0; j < nresonances; ++j){
				tass << "\t\t" << scientific << setprecision(7) << energy[i][j] << defaultfloat << "\t" << gamma0[i][j] << "\t" << gamma[i][j] << "\t" << ji[i] << "\t" << jj[i][j] << "\n";
			}
			
			// Reset precision
			tass << defaultfloat << setprecision((int) default_precision);
		} else{
			cout << "Error: " << __FILE__ << ":" << __LINE__ << ": " << endl; 
			cout << " printTarget(): Number of parameters ENERGY, GAMMA0, GAMMA, J0 and J does not match."  << endl;
			abort();
		}
	}

	tass << "DOPPLER_BROADENING:\t";

	if(dopplerBroadening.size() == 0){
		tass << "not set" ;
	} else{
		switch(dopplerBroadening[i]){
			case dopplerModel::zero:
				tass << "zero" ;
				break;
			case dopplerModel::arb_vdist:
				tass << "arb_vdist, bins " << velocityBinFile[i] << ", velocity distribution = " << vDistFile[i] ;
				break;
			case dopplerModel::arb_cs:
				tass << "arb_cs, bins = " << energyBinFile[i] << ", cross section = " << crosssectionFile[i] ;
				break;
			case dopplerModel::mb:
				tass << "Maxwell-Boltzmann, T_eff = " << dopplerParams[i][0] << " K" ;
				break;
			case dopplerModel::mba:
				tass << "Maxwell-Boltzmann (using approximation), T_eff = " << dopplerParams[i][0] << " K " ;
				break;
			case dopplerModel::mbd:
				tass << "Maxwell-Boltzmann (using Debye approximation), T = " << dopplerParams[i][0] << " K, T_D = " << dopplerParams[i][1] << " K " ;
				break;
			case dopplerModel::mbad:
				tass << "Maxwell-Boltzmann (using Debye and integral approximation), T = " << dopplerParams[i][0] << " K, T_D = " << dopplerParams[i][1] << " K " ;
				break;
			case dopplerModel::phdos:
				//tass << "phDOS, omega_s = " << omegaFile[i] << ", e_s = " << polarizationFile[i] << ", p = " << momentumFile[i] << ", T = " << dopplerParams[i][0] << " K, N = " << dopplerParams[i][1] ;
				tass << "phDOS, omega_s = " << omegaFile[i] << ", T = " << dopplerParams[i][0] << " K" ;
				break;
			default: break;
		}
	}
	
	tass << "\nATOMIC MASS:\t\t";
	if(mass.size() == 0){
		tass << "not set" ;
	} else{
		tass << mass[i] << " u" ;
	}
	
	tass << "\nMASS ATTENUATION:\t";

	if(mAttParams.size() == 0 && mAttFile.size() == 0){
		tass << "not set" ;
	} else{
		switch(mAtt[i]){
			case mAttModel::constant:
				tass << "constant, " << mAttParams[i][0] ;
				break;
			case mAttModel::nist:
				tass << "nist, " << mAttFile[i] ;
				break;
			case mAttModel::arb:
				tass << "arb" << mAttFile[i] ;
				break;
			default: break;
		}
	}
	
	tass << "\nTARGET THICKNESS:\t";
	if(thickness.size() == 0){
		tass << "not set" ;
	} else{
	 	tass << thickness[i] << " atoms/fm^2" ;
	}

	tass << "\nVELOCITY:\t\t";
	if(velocity.size() == 0){
		tass << "not set" ;
	} else{
		tass << velocity[i] << " m/s" ;
	}
	tass << "\n";
	
	return tass.str();
}

void Settings::printOptions(){
	cout << option_string();
}

void Settings::printExperiment(){
	cout << experiment_string();
}

void Settings::printTarget(unsigned int i){
	cout << target_string(i);
}

void Settings::write_output() const {

	writeOptions();
	writeExperiment();
	unsigned int ntargets = (unsigned int) targetNames.size();
	for(unsigned int i = 0; i < ntargets; ++i)
		writeTarget(i);

}

void Settings::writeOptions() const{

	ofstream ofile;
	ofile.open(outputfile, std::ios_base::out);

        if(!ofile.is_open()){
                cout << "Error: " << __FILE__ << ":" << __LINE__ << ": "; 
		cout << " writeOptions(): File '" << outputfile << "' could not be opened." << endl;
		abort();
	}

	ofile << option_string();
	ofile.close();
}

void Settings::writeExperiment() const{

	ofstream ofile;
	ofile.open(outputfile, std::ios_base::out | std::ios_base::app);

        if(!ofile.is_open()){
                cout << "Error: " << __FILE__ << ":" << __LINE__ << ": "; 
		cout << " write_input(): File '" << outputfile << "' could not be opened." << endl;
		abort();
	}

	ofile << experiment_string();
	ofile.close();
}

void Settings::writeTarget(unsigned int i) const{

	ofstream ofile;
	ofile.open(outputfile, std::ios_base::out | std::ios_base::app);

        if(!ofile.is_open()){
                cout << "Error: " << __FILE__ << ":" << __LINE__ << ": "; 
		cout << " write_input(): File '" << outputfile << "' could not be opened." << endl;
		abort();
	}

	ofile << target_string(i);
	ofile.close();
}
