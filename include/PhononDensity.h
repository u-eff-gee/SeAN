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

#include "Config.h"
#include "Settings.h"
#include "Integrator.h"
#include "Writer.h"

using std::vector;

class PhononDensity{

private:
	vector<double> mu_bins;
	vector<double> mu_hist;
	vector<double> alpha_s_mean;
	vector<double> sine_sum;
	vector<double> cosine_sum;
	vector<double> q_s_squared_over_E_squared;

	Settings settings;
	Integrator integrator;
	Writer writer;

public:
	PhononDensity(Settings &s): integrator(), writer(s){ 
		settings = s; 
	};

	~PhononDensity(){}

	void calculate_mu_bins(const unsigned int target_number);
	void calculate_alpha_s_mean(const vector<double> &omega_s_file, const unsigned int target_number, const unsigned int n_modes);
	void calculate_q_s_squared_over_E_squared(const vector<double> &omega_s_file, const unsigned int target_number, const unsigned int n_modes);

	void calculate_sine_sum(const vector<double> &omega_s_file, const unsigned int n_modes);
	void calculate_cosine_sum(const vector<double> &omega_s_file, const unsigned int n_modes);

	void calculateCrossSection(const vector<double> &energy_bins, const vector<double> &energy_boosted, vector<double> &crosssection_histogram, const vector<double> &omega_s_file, const unsigned int target_number);

// **********************************************************************
// Idealistic implementation, if ALL phonons were known
// **********************************************************************
//
//	void calculate_mu_bins(const unsigned int target_number);
//	void calculate_alpha_s_mean(const vector<double> &omega_s_file, const unsigned int target_number, const unsigned int n_modes);
//	void calculate_q_s_squared_over_E_squared(const vector<double> &omega_s_file, const vector< vector<double> > &e_s_file, const double p[3], const unsigned int target_number, const unsigned int n_modes);
//
//	void calculate_sine_sum(const vector<double> &omega_s_file, const unsigned int n_modes);
//	void calculate_cosine_sum(const vector<double> &omega_s_file, const unsigned int n_modes);
//
//	void calculateCrossSection(const vector<double> &energy_bins, const vector<double> &energy_boosted, vector<double> &crosssection_histogram, const vector<double> &omega_s_file, const vector< vector<double> > (&e_s_file), const vector< vector<double> > (&p_file), const unsigned int target_number);
};

#endif
