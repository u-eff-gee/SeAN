#include "PhononDensity.h"

#include "omp.h"

void PhononDensity::calculateCrossSection(const vector<double> &energy_bins, vector<double> &crosssection_histogram, const vector<double> &omega_s_file, const vector< vector<double> > &e_s_file, const vector< vector<double> > &p_file, const unsigned int target_number){
	n_modes = omega_s_file.size();
	calculate_mu_bins();
	calculate_alpha_s_mean(omega_s_file, n_modes);

	double p[3] = {0., 0., 0.};

	mu_hist = vector<double>(settings.nbins_e);

	for(unsigned int i = 0; i < p_file[0].size(), ++i){
		p[0] = p_file[i][0];
		p[1] = p_file[i][1];
		p[2] = p_file[i][2];

		calculate_q_s_squared_over_e_squared(omega_s_file, e_s_file, p, n_modes);
		calculate_sine_sum(omega_s_file, n_modes);
		calculate_cosine_sum(omega_s_file, n_modes);

		for(unsigned j = 0; j < settings.nbins_e; ++j){
			#pragma omp parallel for
			for(unsigned int k = 0; k < settings.nbins_e; ++k){
				mu_hist[k] = cos(mu_bins[k]*(energy_bins[j] - settings.energy[target_number][0]) - energy_bins[j]*energy_bins[j]*sine_sum[k])*exp(-mu_bins[k]*0.5*settings.gamma[target_number][0] + energy_bins[j]*energy_bins[j]*cosine_sum[k]);
			}
			
			crosssection_histogram[j] += integrator.trapezoidal_rule(mu_bins, mu_hist);
		}
	}

	double cs_max = PI*0.5*HBARC2/(energy_boosted[i]*energy_boosted[i])*(2.*settings.jj[target_number][0] + 1.)/(2.*settings.ji[target_number] + 1.)*settings.gamma0[target_number][0];

	#pragma omp parallel for
	for(unsigned int i = 0; i < settings.nbins_e; ++i){
		crosssection_histogram[i] = cs_max*crosssection_histogram[i];
	}
}

void PhononDensity::calculate_mu_bins(){

	// Calculate integration range for mu. The integrated function is proportional to exp(-mu*Gamma/2). This term will eventually dominate the behavior for large mu, so integrate over a range from 0 to 3 times the "decay constant" 2/Gamma.
	double decay_constant = 2./settings.gamma[target_number][0];
	mu_bins = vector<double>(settings.nbins_e);
	bin_size_mu = (double) 3.*decay_constant/(settings.nbins_e - 1);

	#pragma omp parallel for
	for unsigned i = 0; i < settings.nbins_e; ++i){
		mu_bins[i] = i*bin_size_mu;
	}
}

void PhononDensity::calculate_alpha_s_mean(const vector<double> &omega_s_file, const unsigned int n_modes){
	alpha_s_mean = vector<double>(n_modes);

	double inverse_kBT = 1./(kB * settings.dopplerParams[target_number][0]);

	#pragma omp parallel for
	for(unsigned int i = 0; i < n_modes; ++i){
		alpha_s_mean[i] = 1./(exp(-omega_s_file[i]*inverse_kBT) - 1.);
	}
}

void PhononDensity::calculate_q_s_squared_over_e_squared(const vector<double> &omega_s_file, const vector< vector<double> > &e_s_file, const double &p, const unsigned int n_modes){

	q_s_squared_over_E_squared = vector<double>(n_modes, 0.);

	double c1 = 1./(2.*settings.mass[target_number]*n_modes);

	for(unsigned int i = 0; i < n_modes; ++i){
		q_s_squared_over_E_squared[i] = pow(p[0]*e_s_file[i][0] + p[1]*e_s_file[i][1] + p[2]*e_s_file[i][2], 2)*c1/omega_s_file[i];
	}
}

void PhononDensity::calculate_sine_sum(const vector<double> &omega_s_file, const unsigned int n_modes){

	sine_sum = vector<double>(settings.nbins_e, 0.);

	#pragma omp parallel for
	for(unsigned int i = 0; i < settings.nbins_e; ++i){
		for(unsigned int j = 0; j < n_modes; ++j){
			sine_sum[i] = q_s_squared_over_E_squared[j]*sin(mu_bins[i]*omega_s_file[j]);
		}
	}
}

void PhononDensity::calculate_cosine_sum(const vector<double> &omega_s_file, const unsigned int n_modes){

	cosine_sum = vector<double>(settings.nbins_e, 0.);

	#pragma omp parallel for
	for(unsigned int i = 0; i < settings.nbins_e; ++i){
		for(unsigned int j = 0; j < n_modes; ++j){
			cosine_sum[i] = q_s_squared_over_E_squared[j]*(2.*alpha_s_mean[j] + 1.)*cos(mu_bins[i]*omega_s_file[j]);
		}
	}
}
