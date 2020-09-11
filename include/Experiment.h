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

#include <vector>
#include <iostream>

#include "Target.h"
#include "Settings.h"

using std::vector;
using std::endl;
using std::cout;

class Experiment{

private:
	vector<double> energy_bins;
	vector<Target> targets;

	Settings settings;
	Writer writer;

public:
	//Experiment(){};
	Experiment(Settings &s);
	~Experiment() = default;

	// Functions to manage the calculation process
	void initialize();
	void crossSections();
	void transmission();
	void transmission_thin_target();
	void resonant_scattering();
	void resonant_scattering_thin_target();

	// Functions for output
	string result_string(unsigned int n_setting) const;
	string uncertainty_string() const;
	void plot(unsigned int n_setting) ;
	void print_results(unsigned int n_setting);
	void write(unsigned int n_setting) ;
	void write_results(string outputfile, unsigned int n_setting) const;

private:
	void createEnergyBins(double emin, double emax);
	void createTargets();
};