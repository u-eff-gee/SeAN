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


#include "PhononDensity.h"

#include "omp.h"
#include <cmath>
#include <iostream>
#include <sstream>

using std::cout;
using std::endl;
using std::stringstream;

// **********************************************************************
// Approximation p*e_s ~ 1/3
// An implementation that uses the true eigenvectors e_s and averages
// over a set of vectors p can be found below
// **********************************************************************


void PhononDensity::calculateCrossSection(const vector<double> &energy_bins, const vector<double> &energy_boosted, vector<double> &crosssection_histogram, const vector<double> &omega_s_file, const unsigned int target_number){

	unsigned int n_modes = (unsigned int) omega_s_file.size();

	if(settings.verbosity > 0){
		cout << HORIZONTAL_LINE << endl;
		cout << "> Calculating cross section from " << n_modes << " phonon modes of crystal lattice." << endl; 
	}

	calculate_mu_bins(target_number);

	calculate_alpha_s_mean(omega_s_file, target_number, n_modes);

	mu_hist = vector<double>(settings.nbins_e);

	stringstream save_progress_filename;

	calculate_q_s_squared_over_E_squared(omega_s_file, target_number, n_modes);
	calculate_sine_sum(omega_s_file, n_modes);
	calculate_cosine_sum(omega_s_file, n_modes);

	for(unsigned j = 0; j < settings.nbins_e; ++j){
		#pragma omp parallel for
		for(unsigned int k = 0; k < settings.nbins_e; ++k){
			mu_hist[k] = cos(mu_bins[k]*(energy_bins[j] - energy_boosted[0]) - energy_bins[j]*energy_bins[j]*sine_sum[k])*exp(-mu_bins[k]*0.5*settings.gamma[target_number][0] + energy_bins[j]*energy_bins[j]*cosine_sum[k]);
		}
		
		// Trapezoidal integral seems to be absolutely necessary for this highly oscillating function
		crosssection_histogram[j] += integrator.trapezoidal_rule(mu_bins, mu_hist);

// Debugging information
//		if(settings.verbosity > 1){
//			cout << "\t> Checking normalization : integral W(E) dE = " << integrator.trapezoidal_rule(energy_bins, crosssection_histogram) << " == " << (0.5*PI - atan(-2.*energy_boosted[0]/settings.gamma[target_number][0])) << " ? " << endl;
//		}
	}
	

	double cs_max = PI*HBARC2/(energy_boosted[0]*energy_boosted[0])*(2.*settings.jj[target_number][0] + 1.)/(2.*settings.ji[target_number] + 1.)*settings.gamma0[target_number][0];

	#pragma omp parallel for
	for(unsigned int i = 0; i < settings.nbins_e; ++i){
		crosssection_histogram[i] = cs_max*crosssection_histogram[i];
	}
}

void PhononDensity::calculate_mu_bins(const unsigned int target_number){

	// Calculate integration range for mu. The integrated function is proportional to exp(-mu*Gamma/2). This term will eventually dominate the behavior for large mu, so integrate over a range from 0 to MU_MAX_INTEGRAL times the "decay constant" 2/Gamma.
	double decay_constant = 2./settings.gamma[target_number][0];
	mu_bins = vector<double>(settings.nbins_e);
	double bin_size_mu = (double) MU_MAX_INTEGRAL*decay_constant/(settings.nbins_e - 1);

	#pragma omp parallel for
	for(unsigned i = 0; i < settings.nbins_e; ++i){
		mu_bins[i] = i*bin_size_mu;
	}
}

void PhononDensity::calculate_alpha_s_mean(const vector<double> &omega_s_file, const unsigned int target_number, const unsigned int n_modes){
	alpha_s_mean = vector<double>(n_modes);

	double inverse_kBT = 1./(kB * settings.dopplerParams[target_number][0]);

	#pragma omp parallel for
	for(unsigned int i = 0; i < n_modes; ++i){
		alpha_s_mean[i] = 1./(exp(omega_s_file[i]*inverse_kBT) - 1.);
	}
}

void PhononDensity::calculate_q_s_squared_over_E_squared(const vector<double> &omega_s_file, const unsigned int target_number, const unsigned int n_modes){

	q_s_squared_over_E_squared = vector<double>(n_modes, 0.);

	// Use the fact that the mean value of p*e_s = 1/3
	double c1 = 1./(2.*settings.mass[target_number]*AtomicMassUnit*n_modes);

	#pragma omp parallel for
	for(unsigned int i = 0; i < n_modes; ++i){
		q_s_squared_over_E_squared[i] = c1/omega_s_file[i];
	}
}

void PhononDensity::calculate_sine_sum(const vector<double> &omega_s_file, const unsigned int n_modes){

	sine_sum = vector<double>(settings.nbins_e, 0.);

	#pragma omp parallel for
	for(unsigned int i = 0; i < settings.nbins_e; ++i){
		for(unsigned int j = 0; j < n_modes; ++j){
			sine_sum[i] += q_s_squared_over_E_squared[j]*sin(mu_bins[i]*omega_s_file[j]);
		}
	}
}

void PhononDensity::calculate_cosine_sum(const vector<double> &omega_s_file, const unsigned int n_modes){

	cosine_sum = vector<double>(settings.nbins_e, 0.);

	#pragma omp parallel for
	for(unsigned int i = 0; i < settings.nbins_e; ++i){
		for(unsigned int j = 0; j < n_modes; ++j){
			cosine_sum[i] += q_s_squared_over_E_squared[j]*(2.*alpha_s_mean[j] + 1.)*(cos(mu_bins[i]*omega_s_file[j]) - 1.);
		}
	}
}

// **********************************************************************
// Idealistic implementation, if ALL phonons were known
// **********************************************************************
//
//void PhononDensity::calculateCrossSection(const vector<double> &energy_bins, const vector<double> &energy_boosted, vector<double> &crosssection_histogram, const vector<double> &omega_s_file, const vector< vector<double> > &e_s_file, const vector< vector<double> > &p_file, const unsigned int target_number){
//
//	unsigned int n_modes = (unsigned int) omega_s_file.size();
//	unsigned int n_momentum_vectors = (unsigned int) p_file[0].size();
//	double momentum_vector_weight = (double) 1./n_momentum_vectors;
//
//	if(settings.verbosity > 0){
//		cout << HORIZONTAL_LINE << endl;
//		cout << "> Calculating cross section from " << n_modes << " phonon modes of crystal lattice."; 
//		if(p_file[0].size() > 1){
//			cout << " Averaging over " << n_momentum_vectors << " momentum vectors of the photon." << endl;
//		} else{
//			cout << endl;
//		}
//	}
//
//	calculate_mu_bins(target_number);
//
//	calculate_alpha_s_mean(omega_s_file, target_number, n_modes);
//
//	double p[3] = {0., 0., 0.};
//
//	mu_hist = vector<double>(settings.nbins_e);
//
//	stringstream save_progress_filename;
//
//	for(unsigned int i = 0; i < n_momentum_vectors; ++i){
//		p[0] = p_file[0][i];
//		p[1] = p_file[1][i];
//		p[2] = p_file[2][i];
//
//		if(settings.verbosity > 0){
//			cout << "\tMomentum vector #" << i << " of " << n_momentum_vectors << " : ( " << p[0] << ", " << p[1] << ", " << p[2] << " )" << endl;
//		}
//
//		calculate_q_s_squared_over_E_squared(omega_s_file, e_s_file, p, target_number, n_modes);
//		calculate_sine_sum(omega_s_file, n_modes);
//		calculate_cosine_sum(omega_s_file, n_modes);
//
//		for(unsigned j = 0; j < settings.nbins_e; ++j){
//			#pragma omp parallel for
//			for(unsigned int k = 0; k < settings.nbins_e; ++k){
//				mu_hist[k] = cos(mu_bins[k]*(energy_bins[j] - energy_boosted[0]) - energy_bins[j]*energy_bins[j]*sine_sum[k])*exp(-mu_bins[k]*0.5*settings.gamma[target_number][0] + energy_bins[j]*energy_bins[j]*cosine_sum[k]);
//			}
//			
//			// Trapezoidal integral seems to be absolutely necessary for this highly oscillating function
//			crosssection_histogram[j] += integrator.trapezoidal_rule(mu_bins, mu_hist);
//		}
//
//		if(settings.verbosity > 0){
//			cout << "\t> Checking normalization : integral W(E) dE = " << integrator.trapezoidal_rule(energy_bins, crosssection_histogram) << " == " << (i+1)*(0.5*PI - atan(-2.*energy_boosted[0]/settings.gamma[target_number][0])) << " ? " << endl;
//		}
//
//		if(i % SAVE_PROGRESS == 0){
//			save_progress_filename << settings.targetNames[0] << "_crosssection_save_" << i;
//			writer.write1DHistogram(crosssection_histogram, save_progress_filename.str(), "Cross section / fm^2");
//			save_progress_filename.str("");
//			save_progress_filename.clear();
//		}
//	}
//	
//
//	double cs_max = PI*HBARC2/(energy_boosted[0]*energy_boosted[0])*(2.*settings.jj[target_number][0] + 1.)/(2.*settings.ji[target_number] + 1.)*settings.gamma0[target_number][0]*momentum_vector_weight;
//
//	#pragma omp parallel for
//	for(unsigned int i = 0; i < settings.nbins_e; ++i){
//		crosssection_histogram[i] = cs_max*crosssection_histogram[i];
//	}
//}
//
//void PhononDensity::calculate_mu_bins(const unsigned int target_number){
//
//	// Calculate integration range for mu. The integrated function is proportional to exp(-mu*Gamma/2). This term will eventually dominate the behavior for large mu, so integrate over a range from 0 to MU_MAX_INTEGRAL times the "decay constant" 2/Gamma.
//	double decay_constant = 2./settings.gamma[target_number][0];
//	mu_bins = vector<double>(settings.nbins_e);
//	double bin_size_mu = (double) MU_MAX_INTEGRAL*decay_constant/(settings.nbins_e - 1);
//
//	#pragma omp parallel for
//	for(unsigned i = 0; i < settings.nbins_e; ++i){
//		mu_bins[i] = i*bin_size_mu;
//	}
//}
//
//void PhononDensity::calculate_alpha_s_mean(const vector<double> &omega_s_file, const unsigned int target_number, const unsigned int n_modes){
//	alpha_s_mean = vector<double>(n_modes);
//
//	double inverse_kBT = 1./(kB * settings.dopplerParams[target_number][0]);
//
//	#pragma omp parallel for
//	for(unsigned int i = 0; i < n_modes; ++i){
//		alpha_s_mean[i] = 1./(exp(omega_s_file[i]*inverse_kBT) - 1.);
//	}
//}
//
//void PhononDensity::calculate_q_s_squared_over_E_squared(const vector<double> &omega_s_file, const vector< vector<double> > &e_s_file, const double p[3], const unsigned int target_number, const unsigned int n_modes){
//
//	q_s_squared_over_E_squared = vector<double>(n_modes, 0.);
//
//	//double c1 = 1./(2.*settings.mass[target_number]*AtomicMassUnit*n_modes);
//	// Use the fact that the mean value of p*e_s = 1/3
//	double c1 = 1./(2.*settings.mass[target_number]*AtomicMassUnit*n_modes);
//
//	#pragma omp parallel for
//	for(unsigned int i = 0; i < n_modes; ++i){
//		//q_s_squared_over_E_squared[i] = pow(p[0]*e_s_file[0][i] + p[1]*e_s_file[1][i] + p[2]*e_s_file[2][i], 2)*c1/omega_s_file[i];
//		q_s_squared_over_E_squared[i] = c1/omega_s_file[i];
//	}
//}
//
//void PhononDensity::calculate_sine_sum(const vector<double> &omega_s_file, const unsigned int n_modes){
//
//	sine_sum = vector<double>(settings.nbins_e, 0.);
//
//	#pragma omp parallel for
//	for(unsigned int i = 0; i < settings.nbins_e; ++i){
//		for(unsigned int j = 0; j < n_modes; ++j){
//			sine_sum[i] += q_s_squared_over_E_squared[j]*sin(mu_bins[i]*omega_s_file[j]);
//		}
//	}
//}
//
//void PhononDensity::calculate_cosine_sum(const vector<double> &omega_s_file, const unsigned int n_modes){
//
//	cosine_sum = vector<double>(settings.nbins_e, 0.);
//
//	#pragma omp parallel for
//	for(unsigned int i = 0; i < settings.nbins_e; ++i){
//		for(unsigned int j = 0; j < n_modes; ++j){
//			cosine_sum[i] += q_s_squared_over_E_squared[j]*(2.*alpha_s_mean[j] + 1.)*(cos(mu_bins[i]*omega_s_file[j]) - 1.);
//		}
//	}
//}
