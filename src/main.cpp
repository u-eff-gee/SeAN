#include <iostream>
#include <argp.h>
#include <ctime>
#include <chrono>

#include "Experiment.h"
#include "Config.h"
#include "Settings.h"
#include "InputReader.h"

using std::cout;
using std::endl;
using std::cin;
using namespace std::chrono;

//const char *sean_program_version = "SeAN 0.0.0";
//const char *sean_program_bug_address = "<ugayer@ikp.tu-darmstadt.de>";

static char doc[] = "SeAN, Self-Absorption Numerical";
static char args_doc[] = "INPUTFILE";

static struct argp_option options[] = {
  { 0, 'e', 0, 0, "Do not use convolution approximation (increased computing time)" },
  { 0, 'p', 0, 0, "Create plots of all calculated quantities" },
  { 0, 'w', 0, 0, "Create text output files for all calculated quantities" },
  { 0 }
};

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
    struct Settings *settings = (struct Settings*) state->input;

    switch (key) {
    	case ARGP_KEY_ARG: settings->inputfile = arg; break;
	case 'e': settings->exact = true; break;
    	case 'p': settings->plot = true; break;
    	case 'w': settings->write = true; break;
    	case ARGP_KEY_END:
        	if(state->arg_num == 0) {
            		argp_usage(state);
        	}
        	if(settings->inputfile == "") {
            		argp_failure(state, 1, 0, "INPUTFILE not specified.");
        	}
        	break;
    	default: return ARGP_ERR_UNKNOWN;
    	}
    	return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc, 0, 0, 0 };

int main(int argc, char* argv[]){

	vector <Settings> settings;
	settings.push_back(Settings());
	argp_parse(&argp, argc, argv, 0, 0, &settings[0]);

	// Start the clock
	high_resolution_clock::time_point start = high_resolution_clock::now();


	if(settings[0].sudowrite)
		settings[0].write = true;

	InputReader input;
	input.readFile(settings);
	
	for(auto s: settings)
		s.print();

//	Experiment experiment(settings);
//	experiment.initialize();
//	experiment.crossSections();
//	experiment.transmission();
//	if(settings.plot){
//		experiment.plot();
//	}
//	if(settings.write){
//		experiment.write();
//	}
//	experiment.resonant_scattering();
//	experiment.print_results();
//	
//	// Stop the clock
//	high_resolution_clock::time_point stop = high_resolution_clock::now();
//	duration<double> delta_t = duration_cast< duration<double>>(stop - start);
//
//	cout << "> main.cpp: Execution took " << delta_t.count() << " seconds" << endl;
}
