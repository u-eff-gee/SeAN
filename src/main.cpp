#include <iostream>
#include <argp.h>
#include <ctime>
#include <chrono>
#include <sstream>

#include "Experiment.h"
#include "Config.h"
#include "Settings.h"
#include "InputReader.h"

using std::cout;
using std::endl;
using std::cin;
using std::stringstream;
using namespace std::chrono;

//const char *sean_program_version = "SeAN 0.0.0";
//const char *sean_program_bug_address = "<ugayer@ikp.tu-darmstadt.de>";

static char doc[] = "SeAN, Self-Absorption Numerical";
static char args_doc[] = "INPUTFILE";

static struct argp_option options[] = {
  { "exact", 'e', 0, 0, "Do not use convolution approximation (increased computing time)" },
  { "plot", 'p', 0, 0, "Create plots of all calculated quantities" },
  { "write", 'w', 0, 0, "Create text output files for all calculated quantities" },
  { "verbosity", 'v', "VERBOSITY", 0, "Set command line verbosity (0 = print nothing, 1 = print results, 2 [default] = print input and results)" },
  { "recoil", 'r', 0, 0, "Include nuclear recoil in the calculation of the cross section maximum."},
  { "output", 'o', "OUTPUTFILENAME", 0, "Write input and results to file" },
  { 0 }
};

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
    struct Settings *settings = (struct Settings*) state->input;

    switch (key) {
    	case ARGP_KEY_ARG: settings->inputfile = arg; break;
	case 'e': settings->exact = true; break;
    	case 'p': settings->plot = true; break;
    	case 'w': settings->write = true; break;
	case 'v': settings->verbosity = atoi(arg); break;
	case 'r': settings->recoil= true; break;
	case 'o': settings->output = true;
		  settings->outputfile = arg;
			break;
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

	InputReader input;
	input.readFile(settings);

	Experiment *experiment;

	long unsigned int n_settings = settings.size();

	for(unsigned int i = 0; i < n_settings; ++i){

		if(settings[0].verbosity > 1){
			settings[i].print();
		}
		experiment = new Experiment(settings[i]);
		experiment->initialize();
		experiment->crossSections();
		experiment->transmission();
		if(settings[0].plot){
			experiment->plot(i);
		}
		if(settings[0].write){
			experiment->write(i);
		}
		experiment->resonant_scattering();
		if(settings[0].verbosity > 0){
			experiment->print_results();
		}

		if(settings[0].output){
			settings[i].write_output(i);
			experiment->write_results(settings[0].outputfile, i);
		}
	}
	
	if(settings[0].verbosity > 0 && settings[0].output){
		cout << "> Created output file 'output/" << settings[0].outputfile << "'" << endl;
	}
	
	// Stop the clock
	high_resolution_clock::time_point stop = high_resolution_clock::now();
	duration<double> delta_t = duration_cast< duration<double>>(stop - start);

	if(settings[0].verbosity > 0){
		cout << "> main.cpp: Execution took " << delta_t.count() << " seconds" << endl;
	}
}
