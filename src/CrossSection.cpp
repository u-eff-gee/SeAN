#include "CrossSection.h"

void CrossSection::breit_wigner(double (&energies)[NPARAMETERS], double (&crosssection)[NPARAMETERS], vector<double> &e0, vector<double> &gamma0, vector<double> &gamma, vector<double> &jj, double j0){

	double cs_max = 0.;

	for(unsigned int i = 0; i < e0.size(); ++i){
		cs_max = PI*0.5*HBARC2/(e0[i]*e0[i])*(2.*jj[i] + 1.)/(2. * j0 + 1.)*gamma0[i]*gamma[i];
		for(int j = 0; j < NBINS; ++j){
			crosssection[i] += cs_max / ((energies[j] - e0[i])*(energies[j] - e0[i]) + 0.25*gamma[i]*gamma[i]);
		}
	}
}

//
//void CrossSection::maxwellAverage(double (&energies)[NPARAMETERS], double (&crosssection)[NPARAMETERS]){
//	;	
//}
