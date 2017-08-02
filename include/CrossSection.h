#ifndef CROSSSECTION_H 
#define CROSSSECTION_H 1

#include <vector>

#include "Config.h"
#include "math.h"

using std::vector;

class CrossSection{
public:
	CrossSection(){};
	
	~CrossSection(){};

	void breit_wigner(double (&energies)[NPARAMETERS], double (&crosssection)[NPARAMETERS], vector<double> &e0, vector<double> &gamma0, vector<double> &gamma, vector<double> &jj, double j0);

	void absolute_zero(){};

	void maxwell_boltzmann(double (&energies)[NPARAMETERS], double (&crosssection)[NPARAMETERS], vector<double> params);

	void maxwell_boltzmann_debye(double (&energies)[NPARAMETERS], double (&crosssection)[NPARAMETERS], vector<double> params);

private:
	double eGamma(double energy, double velocity){
		return (1. + velocity)/(sqrt(1. - velocity*velocity))*energy;
	}
};

#endif
