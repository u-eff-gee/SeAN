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
	{ "direct", 'd', 0, 0, "Direct calculation of self-absorption. Avoids storing some intermediate quantities, which make SeAN less demanding on memory (default: false).", 0 },
	{ "exact", 'e', 0, 0, "Do not use convolution approximation (increased computing time, default: false)", 0 },
	{ "multi", 'm', 0, 0, "Compute multidimensional integrals (increased computing time, default: false)", 0 },
	{ "output", 'o', "OUTPUTFILENAME", 0, "Write input and results to a file called OUTPUTFILENAME (default: no output writing)", 0 },
	{ "plot", 'p', 0, 0, "Create plots of all calculated quantities (default: false).", 0 },
	{ "recoil", 'r', 0, 0, "Include nuclear recoil in the calculation of the cross section maximum (default: false).", 0 },
	{ "uncertainty", 'u', 0, 0, "Estimate the uncertainty of the numerical evaluations (default: false).", 0 },
	{ "verbosity", 'v', "VERBOSITY", 0, "Set command line verbosity (0 = print nothing, 1 = print results, 2 [default] = print input and results)", 0 },
	{ "write", 'w', 0, 0, "Create text output files for all calculated quantities (default: false).", 0 },
	{ 0, 0, 0, 0, 0, 0 }
};

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
    struct Settings *settings = (struct Settings*) state->input;

    switch (key) {
    case ARGP_KEY_ARG: settings->inputfile = arg; break;
	case 'd': settings->direct = true; break;
	case 'e': settings->exact = true; break;
	case 'm': settings->multi = true; break;
	case 'o': settings->output = true;
		  settings->outputfile = arg;
			break;
    	case 'p': settings->plot = true; break;
	case 'r': settings->recoil= true; break;
	case 'u': settings->uncertainty=true; break;
	case 'v': settings->verbosity = atoi(arg); break;
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

	InputReader input;
	input.readFile(settings);

	long unsigned int n_settings = settings.size();

	for(unsigned int i = 0; i < n_settings; ++i){

		if(settings[0].verbosity > 1){
			settings[i].print();
		}
		Experiment experiment = Experiment(settings[i]);
		experiment.initialize();
		experiment.crossSections();
		experiment.transmission();
		if(settings[0].plot){
			experiment.plot(i);
		}
		if(settings[0].write){
			experiment.write(i);
		}
		experiment.resonant_scattering();
		if(settings[0].verbosity > 0){
			experiment.print_results(i);
		}

		if(settings[0].output){
			settings[i].write_output(i);
			experiment.write_results(settings[0].outputfile, i);
		}
	}
	
	if(settings[0].verbosity > 0 && settings[0].output){
		cout << "> Created output file '" << settings[0].outputfile << "'" << endl;
	}
	
	// Stop the clock
	high_resolution_clock::time_point stop = high_resolution_clock::now();
	duration<double> delta_t = duration_cast< duration<double>>(stop - start);

	if(settings[0].verbosity > 0){
		cout << "> main.cpp: Execution took " << delta_t.count() << " seconds" << endl;
	}
}
