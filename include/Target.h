#ifndef TARGET_H 
#define TARGET_H 1

#include <string>
#include <vector>

#include "Config.h"
#include "CrossSection.h"

using std::string;
using std::vector;

class Target{
private:
	vector<double> e0;
	vector<double> gamma0;
	vector<double> gamma;
	vector<double> jj;
	vector<double> vDistParams;

	CrossSection *crossSection;

	double incident_beam[NBINS] = {0.};
	double crosssection[NBINS] = {0.};
	double massattenuation[NBINS] = {0.};
	double transmitted_beam[NBINS] = {0.};

	string target_name;
	string vDist_ID;
	string massAttenuation_ID;
	double j0;
	double mass;
	double z;

	int target_number;

public:	
	Target(string name, int number){
		target_name = name;
		target_number = number;

		crossSection = new CrossSection();
	};
	
	~Target(){
		delete crossSection;
	};

	void addEnergy(double e){ e0.push_back(e); };
	void addGamma0(double g0){ gamma0.push_back(g0); };
	void addGamma(double g){ gamma.push_back(g); };
	void addJJ(double j){ jj.push_back(j); };
	void setJ0(double j){ j0 = j; };
	void setMass(double m){ mass = m; };
	void setZ(double zz){ z = zz; };

	void setVDist(string vDist){vDist_ID = vDist; };
	void addVDistParameter(double p){ vDistParams.push_back(p); };

	void setMassAttenuation(string mAtt){ massAttenuation_ID = mAtt; };

	void print();
};

#endif
