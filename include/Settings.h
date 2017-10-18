#ifndef SETTINGS_H
#define SETTINGS_H 1

struct Settings{
        char *inputfile;
	bool exact = false;
	bool plot = false;
	bool write = false;
	bool sudowrite = false;

	double emin;
	double emax;
	unsigned int nbins_e;
	unsigned int nbins_z;
};

#endif
