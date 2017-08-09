#include <iostream>
#include <argp.h>
#include <ctime>
#include <chrono>

#include "Experiment.h"
#include "Config.h"

using std::cout;
using std::endl;
using namespace std::chrono;

//const char *sean_program_version = "SeAN 0.0.0";
//const char *sean_program_bug_address = "<ugayer@ikp.tu-darmstadt.de>";

static char doc[] = "SeAN, Self-Absorption Numerical";
static char args_doc[] = "INPUTFILE";

static struct argp_option options[] = {
  { 0, 'p', 0, 0, "Create plots of all calculated quantities" },
  { 0, 'o', 0, 0, "Create output files for all calculated quantities" },
  { 0 }
};

struct arguments{
        char *inputfile;
	bool plot = false;
};

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
    struct arguments *arguments = (struct arguments*)state->input;

    switch (key) {
    case ARGP_KEY_ARG: arguments->inputfile = arg; break;
    case 'p': arguments->plot = true; break;
    case ARGP_KEY_END:
        if(state->arg_num == 0) {
            argp_usage(state);
        }
        if(!arguments->inputfile) {
            argp_failure(state, 1, 0, "INPUTFILE not specified.");
        }
        break;
    default: return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc, 0, 0, 0 };

int main(int argc, char* argv[]){

	struct arguments arguments;
	argp_parse(&argp, argc, argv, 0, 0, &arguments);

	high_resolution_clock::time_point start = high_resolution_clock::now();

	Experiment experiment;

	experiment.readInputFile(arguments.inputfile);
	experiment.print();
	experiment.crossSections(arguments.plot);
	experiment.transmission(arguments.plot);

	high_resolution_clock::time_point stop = high_resolution_clock::now();
	duration<double> delta_t = duration_cast< duration<double>>(stop - start);
	cout << "> main.cpp: Execution took " << delta_t.count() << " seconds" << endl;
}
